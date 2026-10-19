import test from "node:test";
import assert from "node:assert/strict";
import { parseAbif, parseFastq, parseInputText, readSequenceFile, qualitySummary, trimRead } from "../sequencing.mjs";
import { alignSequences, reverseComplement, toText } from "../align-core.mjs";

function abifFixture({ sequence = "ACGTAC", callSet = 2, quality = [10, 20, 30, 40, 25, 35], positions = [2, 5, 8, 11, 14, 17] } = {}) {
  const entries = [
    { tag: "PBAS", number: callSet, size: 1, data: Buffer.from(sequence) },
    { tag: "FWO_", number: 1, size: 1, data: Buffer.from("GATC") },
    { tag: "SMPL", number: 1, size: 1, type: 18, data: Buffer.from([4, 116, 101, 115, 116]) }
  ];
  if (quality) entries.push({ tag: "PCON", number: callSet, size: 1, data: Buffer.from(quality) });
  if (positions) {
    const data = Buffer.alloc(positions.length * 2);
    positions.forEach((p, i) => data.writeUInt16BE(p, i * 2));
    entries.push({ tag: "PLOC", number: callSet, size: 2, data });
    for (let i = 0; i < 4; i++) {
      const channel = Buffer.alloc(40);
      for (let j = 0; j < 20; j++) channel.writeInt16BE(j + i * 20, j * 2);
      entries.push({ tag: "DATA", number: i + 9, size: 2, data: channel });
    }
  }
  const dataSize = entries.reduce((sum, e) => sum + (e.data.length > 4 ? e.data.length : 0), 0);
  const buffer = Buffer.alloc(128 + 28 * entries.length + dataSize);
  buffer.write("ABIF"); buffer.writeUInt16BE(101, 4);
  buffer.writeUInt32BE(entries.length, 18); buffer.writeUInt32BE(128, 26);
  let offset = 128 + 28 * entries.length;
  entries.forEach((entry, i) => {
    const at = 128 + i * 28;
    buffer.write(entry.tag, at); buffer.writeUInt32BE(entry.number, at + 4);
    buffer.writeUInt16BE(entry.type || (entry.size === 2 ? 4 : 2), at + 8);
    buffer.writeUInt16BE(entry.size, at + 10);
    buffer.writeUInt32BE(entry.data.length / entry.size, at + 12);
    buffer.writeUInt32BE(entry.data.length, at + 16);
    if (entry.data.length <= 4) entry.data.copy(buffer, at + 20);
    else { buffer.writeUInt32BE(offset, at + 20); entry.data.copy(buffer, offset); offset += entry.data.length; }
  });
  return buffer;
}

test("AB1 reads matching call/quality tags, inline FWO, big-endian peaks and channels", () => {
  const record = parseAbif(abifFixture(), "sample-F.ab1");
  assert.equal(record.sequence, "ACGTAC");
  assert.equal(record.name, "sample-F");
  assert.equal(record.sampleName, "test");
  assert.deepEqual(record.qualities, [10, 20, 30, 40, 25, 35]);
  assert.deepEqual(record.trace.positions, [2, 5, 8, 11, 14, 17]);
  assert.equal(record.trace.channels.G[0], 0);
  assert.equal(record.trace.channels.A[0], 20);
  assert.equal(record.trace.channels.T[0], 40);
  assert.equal(record.trace.channels.C[0], 60);
  assert.equal(qualitySummary(record).q20, 100 * 5 / 6);
});

test("AB1 inline short reads, PBAS1 fallback, and missing quality remain usable", () => {
  const record = parseAbif(abifFixture({ sequence: "ACG", quality: [10, 30, 40], positions: null, callSet: 1 }));
  assert.equal(record.sequence, "ACG");
  assert.deepEqual(record.qualities, [10, 30, 40]);
  assert.equal(record.trace, null);
  assert.equal(parseAbif(abifFixture({ quality: null })).qualities, null);
  const mismatch = parseAbif(abifFixture({ quality: [20] }));
  assert.equal(mismatch.qualities, null);
  assert.match(mismatch.warnings.join(""), /不一致/);
});

test("truncated or corrupt AB1 cannot return fabricated sequence or out-of-range data", () => {
  assert.throws(() => parseAbif(Buffer.from("ABIF")), /有效/);
  const broken = abifFixture();
  broken.writeUInt32BE(0xFFFFFFFF, 26);
  assert.throws(() => parseAbif(broken), /目录不完整/);
  const dataBroken = abifFixture();
  dataBroken.writeUInt32BE(0xFFFFFFFF, 128 + 20);
  assert.throws(() => parseAbif(dataBroken), /PBAS2 数据不完整/);
});

test("FASTQ accepts wrapped reads and quality lines beginning with @ or +", () => {
  const records = parseFastq("@first\nACGT\nAC\n+first\n@@++\nII\n@second\nTGCA\n+\n!5I~\n");
  assert.equal(records.length, 2);
  assert.equal(records[0].sequence, "ACGTAC");
  assert.deepEqual(records[0].qualities, [31, 31, 10, 10, 40, 40]);
  assert.deepEqual(records[1].qualities, [0, 20, 40, 93]);
  assert.throws(() => parseFastq("@a\nACGT\n+\n!!!\n"), /质量字符数/);
  assert.throws(() => parseFastq("@a\nACGT\n+b\nIIII\n"), /名称不一致/);
});

test("SEQ comments, FASTQ file detection and archive rejection give useful results", async () => {
  assert.equal(parseInputText(";ABI export\n>sample\n1 acgt 4\n")[0].sequence, "ACGT");
  const seq = await readSequenceFile(new File(["ACGTN\n"], "sample.seq"));
  assert.equal(seq[0].sequence, "ACGTN");
  const fastq = await readSequenceFile(new File(["@a\nACGT\n+\nIIII\n"], "sample.fastq"));
  assert.deepEqual(fastq[0].qualities, [40, 40, 40, 40]);
  const ab1 = await readSequenceFile(new File([abifFixture()], "sample.ab1"));
  assert.equal(ab1[0].format, "AB1");
  await assert.rejects(() => readSequenceFile(new File(["PKstub"], "result.zip")), /先解压/);
  await assert.rejects(() => readSequenceFile(new File(["image"], "result.pdf")), /原始 AB1/);
});

test("quality trimming is reversible, end-only, and refuses uniformly poor reads", () => {
  const sequence = "NN" + "ACGT".repeat(8) + "NNN";
  const qualities = Array.from(sequence, b => b === "N" ? 3 : 35);
  qualities[17] = 5; // internal low-quality site must remain
  const original = { name: "read", sequence, qualities };
  const trimmed = trimRead(original);
  assert.equal(trimmed.sequence, "ACGT".repeat(8));
  assert.equal(trimmed.qualities[15], 5);
  assert.equal(trimmed.sourceStart, 3);
  assert.equal(trimmed.sourceEnd, 34);
  assert.equal(original.sequence, sequence);
  assert.equal(original.qualities.length, sequence.length);
  assert.throws(() => trimRead({ name: "bad", sequence: "A".repeat(20), qualities: Array(20).fill(5) }), /没有连续/);
  assert.equal(trimRead({ name: "text", sequence: "ACGT" }).sequence, "ACGT");
});

const reference = "TTGACCGTACGATGCCTAGGTCGATCTAGCACTTGACGATGGCACTAGGTACCTGACCTAACGTAGGCTAGCTACGATGCCATCGTACAGGACTCAGT";
const fit = (records, orient = false) => alignSequences(records, orient, () => {}, { mode: "fitting" });

test("Sanger fragments locate inside reference without counting uncovered tails as deletion", async () => {
  const fragment = reference.slice(19, 77);
  const result = await fit([{ name: "ref", sequence: reference }, { name: "read", sequence: fragment }]);
  const stat = result.stats[0];
  assert.equal(stat.refStart, 20); assert.equal(stat.refEnd, 77);
  assert.equal(stat.identity, 100); assert.equal(stat.gaps, 0);
  assert.equal(stat.uncovered, reference.length - fragment.length);
  assert.equal(stat.coveredReference, fragment.length);
  assert.equal(result.rows[1].replace(/-/g, ""), fragment);
  assert.match(toText(result), /参照覆盖 20–77/);
});

test("Sanger reverse reads, internal substitutions and deletions retain reference coordinates", async () => {
  const fragment = reference.slice(19, 77);
  const deleted = fragment.slice(0, 20) + fragment.slice(21);
  const result = await fit([{ name: "ref", sequence: reference }, { name: "reverse", sequence: reverseComplement(deleted) }], true);
  assert.equal(result.orientations[1], "反向互补");
  assert.equal(result.stats[0].gaps, 1);
  assert.equal(result.stats[0].substitutions, 0);
  assert.equal(result.stats[0].variants[0].position, "40");
  assert.equal(result.stats[0].uncovered, reference.length - fragment.length);
});

test("multi-read merged gaps are not false differences for another identical read", async () => {
  const records = [{ name: "ref", sequence: reference },
    { name: "insertion", sequence: reference.slice(0, 40) + "A" + reference.slice(40) },
    { name: "identical", sequence: reference }];
  const result = await alignSequences(records, false);
  assert.equal(result.stats[0].gaps, 1);
  assert.equal(result.stats[1].gaps, 0);
  assert.equal(result.stats[1].identity, 100);
});

test("ambiguous bases cannot outscore a canonical exact fragment", async () => {
  const result = await fit([{ name: "ref", sequence: "NNNNACGTNNNN" }, { name: "read", sequence: "ACGT" }]);
  assert.equal(result.stats[0].refStart, 5);
  assert.equal(result.stats[0].refEnd, 8);
  assert.equal(result.stats[0].matches, 4);
});

test("long references are allowed only in bounded fragment mode", async () => {
  const ref = "A".repeat(3010) + "CGTCAGTCCGATGCTAGT";
  const result = await fit([{ name: "ref", sequence: ref }, { name: "read", sequence: "CGTCAGTCCGATGCTAGT" }]);
  assert.equal(result.stats[0].refStart, 3011);
  await assert.rejects(() => fit([{ name: "ref", sequence: "A".repeat(12001) }, { name: "r", sequence: "ACGT" }]), /12,000/);
});
