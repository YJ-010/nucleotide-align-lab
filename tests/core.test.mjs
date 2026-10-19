import test from "node:test";
import assert from "node:assert/strict";
import { parseSequences, alignSequences, reverseComplement, toAlignedFasta, toClustal, toText } from "../align-core.mjs";

test("FASTA input and a single substitution produce correct coordinates and identity", async () => {
  const records = parseSequences(">参照\nACTGGT\n>样本\nACTGAT\n");
  const result = await alignSequences(records, false);
  assert.deepEqual(result.rows, ["ACTGGT", "ACTGAT"]);
  assert.equal(result.stats[0].substitutions, 1);
  assert.equal(result.stats[0].variants[0].position, "5");
  assert.equal(result.stats[0].identity.toFixed(1), "83.3");
  assert.match(toClustal(result), /^CLUSTAL/);
  assert.match(toText(result), /Seq2 \/ Seq1/);
});

test("ApE style GenBank records parse all ORIGIN lines and stop before features", () => {
  const body = [
    "LOCUS       pExample             18 bp    DNA     circular",
    "FEATURES             Location/Qualifiers",
    "     source          1..18",
    "ORIGIN",
    "        1 atg cct agg",
    "       10 gca tta gcc",
    "//"
  ].join("\n");
  const [record] = parseSequences(body);
  assert.equal(record.name, "pExample");
  assert.equal(record.sequence, "ATGCCTAGGGCATTAGCC");
});

test("CLUSTAL blocks are joined and input gaps are removed for new alignment", () => {
  const body = "CLUSTAL nucleotide alignment\n\nSeq1  ACG-T  4\nSeq2  AC-CT  4\n      **  *\n\nSeq1  GG\nSeq2  GA\n";
  const records = parseSequences(body);
  assert.deepEqual(records.map(r => r.sequence), ["ACGTGG", "ACCTGA"]);
});

test("reverse complement is chosen only when it scores better", async () => {
  const reference = "ATGCTACCGTAA";
  const result = await alignSequences([
    { name: "ref", sequence: reference },
    { name: "minus", sequence: reverseComplement(reference) }
  ], true);
  assert.equal(result.orientations[1], "反向互补");
  assert.equal(result.stats[0].identity, 100);
});

test("center-star merge keeps all rows aligned and preserves ungapped source", async () => {
  const records = parseSequences(">ref\nACGTTACG\n>insertion\nACGCTTACG\n>deletion\nACGTACG\n>another insertion\nACGGTTACG");
  const result = await alignSequences(records, false);
  assert.ok(result.rows.every(row => row.length === result.length));
  result.rows.forEach((row, i) => assert.equal(row.replace(/-/g, ""), records[i].sequence));
  assert.ok(result.stats.some(stat => stat.gaps > 0));
  assert.equal((toAlignedFasta(result).match(/^>/gm) || []).length, 4);
});

test("unsupported characters and oversized requests fail with useful errors", async () => {
  assert.throws(() => parseSequences(">a\nATGZ"), /无效字符/);
  await assert.rejects(() => alignSequences([{ name: "one", sequence: "ATG" }]), /2–8/);
  await assert.rejects(() => alignSequences([
    { name: "a", sequence: "A".repeat(3001) },
    { name: "b", sequence: "A" }
  ]), /3,000/);
});
