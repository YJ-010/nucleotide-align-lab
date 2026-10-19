import { parseSequences } from "./align-core.mjs?v=20261019-clone1";

const DNA = /^[ACGTRYSWKMBDHVN]+$/;
const ascii = bytes => Array.from(bytes, byte => String.fromCharCode(byte)).join("");
const cleanName = name => String(name || "测序样本").replace(/[\r\n>]/g, " ").trim().slice(0, 90);

export function parseFastq(input) {
  const lines = input.replace(/^\uFEFF/, "").replace(/\r\n?/g, "\n").split("\n");
  const records = [];
  let i = 0;
  while (i < lines.length) {
    if (!lines[i].trim()) { i++; continue; }
    if (!lines[i].startsWith("@")) throw new Error(`FASTQ 第 ${i + 1} 行应为 @ 样本名称。`);
    const header = lines[i++].slice(1).trim();
    let sequence = "";
    while (i < lines.length && !lines[i].startsWith("+")) sequence += lines[i++].trim();
    sequence = sequence.toUpperCase().replace(/U/g, "T");
    if (!DNA.test(sequence) || i >= lines.length) throw new Error(`FASTQ“${header}”的序列或 + 分隔行不完整。`);
    const repeatedHeader = lines[i++].slice(1).trim();
    if (repeatedHeader && repeatedHeader !== header) throw new Error(`FASTQ“${header}”两处样本名称不一致。`);
    let quality = "";
    while (i < lines.length && quality.length < sequence.length) quality += lines[i++];
    if (quality.length !== sequence.length || /[^!-~]/.test(quality))
      throw new Error(`FASTQ“${header}”的质量字符数必须等于碱基数，且须为 Phred+33。`);
    records.push({ name: cleanName(header), sequence,
      qualities: Array.from(quality, char => char.charCodeAt(0) - 33), format: "FASTQ", warnings: [] });
    if (records.length > 8) throw new Error("本工具一次最多导入 8 条序列；大批量高通量 FASTQ 请使用专门的分析软件。");
  }
  return records;
}

export function parseInputText(input, fallbackName = "序列 1") {
  const text = input.replace(/^\uFEFF/, "").trim();
  if (!text) return [];
  if (text.startsWith("@")) return parseFastq(text);
  // ABI .seq exports can have semicolon comment lines before their sequence.
  const uncommented = text.split(/\r?\n/).filter(line => !/^\s*[;#]/.test(line)).join("\n");
  return parseSequences(uncommented, fallbackName).map(record => ({ ...record, format: "文本", warnings: [] }));
}

export function parseAbif(input, filename = "测序样本.ab1") {
  const bytes = input instanceof Uint8Array ? input : new Uint8Array(input);
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  if (bytes.length < 34 || ascii(bytes.subarray(0, 4)) !== "ABIF")
    throw new Error(`“${filename}”不是有效的 AB1/ABIF 文件，或文件已损坏。`);
  const version = view.getUint16(4);
  if (version < 100 || version >= 200) throw new Error(`“${filename}”的 ABIF 版本暂不支持。`);
  const count = view.getUint32(18), directory = view.getUint32(26);
  if (!count || count > 10000 || directory < 34 || directory + count * 28 > bytes.length)
    throw new Error(`“${filename}”的 AB1 目录不完整，请重新下载原始文件。`);
  const entries = new Map();
  for (let i = 0; i < count; i++) {
    const at = directory + i * 28;
    const tag = ascii(bytes.subarray(at, at + 4));
    entries.set(tag + view.getUint32(at + 4), {
      type: view.getUint16(at + 8), elementSize: view.getUint16(at + 10),
      count: view.getUint32(at + 12), size: view.getUint32(at + 16), at,
      offset: view.getUint32(at + 20)
    });
  }
  function read(tag) {
    const entry = entries.get(tag);
    if (!entry) return null;
    const offset = entry.size <= 4 ? entry.at + 20 : entry.offset;
    if (offset + entry.size > bytes.length || entry.count * entry.elementSize > entry.size)
      throw new Error(`“${filename}”的 ${tag} 数据不完整。`);
    return { ...entry, data: bytes.subarray(offset, offset + entry.size) };
  }
  function text(tag) {
    const entry = read(tag);
    if (!entry) return "";
    let data = entry.data;
    if (entry.type === 18) data = data.subarray(1, 1 + data[0]);
    return ascii(data).replace(/\0+$/, "");
  }
  function integers(tag) {
    const entry = read(tag);
    if (!entry) return null;
    if (![1, 2, 4].includes(entry.elementSize)) throw new Error(`“${filename}”的 ${tag} 数值格式暂不支持。`);
    const values = new DataView(entry.data.buffer, entry.data.byteOffset, entry.data.byteLength);
    return Array.from({ length: entry.count }, (_, i) => entry.elementSize === 1 ? values.getUint8(i)
      : entry.elementSize === 2 ? values.getUint16(i * 2) : values.getUint32(i * 4));
  }
  // Match sequence, quality and peak-position tags from the SAME base-call set.
  const callSet = entries.has("PBAS2") ? 2 : 1;
  const sequence = text("PBAS" + callSet).toUpperCase();
  if (!sequence) throw new Error(`“${filename}”没有已判读的 DNA 序列；请让测序公司提供含序列的 AB1 或 .seq 文件。`);
  if (!DNA.test(sequence)) throw new Error(`“${filename}”的判读序列含不支持的碱基字符。`);
  const warnings = [];
  let qualities = null, trace = null;
  try {
    qualities = integers("PCON" + callSet);
    if (qualities && qualities.length !== sequence.length) {
      qualities = null; warnings.push("质量值与碱基数不一致，未使用质量值。");
    }
  } catch { warnings.push("质量数据不完整，已保留序列。"); }
  try {
    const order = text("FWO_1");
    const positions = integers("PLOC" + callSet);
    if (positions && positions.length === sequence.length && /^[ACGT]{4}$/.test(order) && new Set(order).size === 4) {
      const channels = {};
      for (let i = 0; i < 4; i++) {
        const entry = read("DATA" + (9 + i));
        if (!entry || entry.elementSize !== 2 || !entry.count) throw new Error("missing trace");
        const values = new DataView(entry.data.buffer, entry.data.byteOffset, entry.data.byteLength);
        channels[order[i]] = Array.from({ length: entry.count }, (_, j) => values.getInt16(j * 2));
      }
      const lengths = Object.values(channels).map(channel => channel.length);
      if (new Set(lengths).size !== 1 || positions.some((p, i) => p >= lengths[0] || (i > 0 && p <= positions[i - 1])))
        throw new Error("invalid peak positions");
      trace = { positions, channels };
    }
  } catch { warnings.push("峰图数据不完整，已保留可读的序列。"); }
  return { name: cleanName(filename.replace(/\.[^.]+$/, "")), sequence, qualities, trace,
    format: "AB1", sourceFile: filename, sampleName: text("SMPL1"), warnings };
}

export async function readSequenceFile(file) {
  if (file.size > 10_000_000) throw new Error(`“${file.name}”超过 10 MB，请选择单条 Sanger 测序文件或提取所需序列。`);
  const bytes = new Uint8Array(await file.arrayBuffer());
  if (ascii(bytes.subarray(0, 4)) === "ABIF" || /\.(ab1|abi)$/i.test(file.name))
    return [parseAbif(bytes, file.name)];
  if ((bytes[0] === 80 && bytes[1] === 75) || (bytes[0] === 31 && bytes[1] === 139))
    throw new Error(`“${file.name}”是压缩包，请先解压，再选择其中的 AB1、SEQ、FASTA 或 FASTQ 文件。`);
  if (/\.(pdf|png|jpg|jpeg|docx?|scf)$/i.test(file.name))
    throw new Error(`“${file.name}”暂不能直接比对，请提供原始 AB1 或导出的序列文本。`);
  if (bytes.length > 2_000_000) throw new Error(`文本文件“${file.name}”超过 2 MB；本工具用于少量基因片段比对。`);
  let body;
  if (bytes[0] === 255 && bytes[1] === 254) body = new TextDecoder("utf-16le").decode(bytes);
  else if (bytes[0] === 254 && bytes[1] === 255) body = new TextDecoder("utf-16be").decode(bytes);
  else {
    try { body = new TextDecoder("utf-8", { fatal: true }).decode(bytes); }
    catch { body = new TextDecoder("gb18030").decode(bytes); }
  }
  if (body.includes("\0")) throw new Error(`“${file.name}”不是可识别的序列文本。`);
  return parseInputText(body, cleanName(file.name.replace(/\.[^.]+$/, "")))
    .map(record => ({ ...record, sourceFile: file.name }));
}

export function qualitySummary(record) {
  const q = record.qualities;
  return { length: record.sequence.length,
    ambiguous: Array.from(record.sequence).filter(base => !/[ACGT]/.test(base)).length,
    q20: q ? 100 * q.filter(value => value >= 20).length / q.length : null };
}

export function trimRead(record, threshold = 20, run = 10) {
  if (!record.qualities || record.qualities.length !== record.sequence.length) return { ...record, trimNote: "无质量值，保留全长" };
  let good = 0, first = -1, end = -1;
  for (let i = 0; i < record.sequence.length; i++) {
    good = record.qualities[i] >= threshold && /[ACGT]/.test(record.sequence[i]) ? good + 1 : 0;
    if (good >= run) {
      if (first < 0) first = i - run + 1;
      end = i + 1;
    }
  }
  if (first < 0) throw new Error(`“${record.name}”没有连续 ${run} 个 Q≥${threshold} 的明确碱基。请检查测序质量，或取消自动修剪后查看原始读段。`);
  return { ...record, sequence: record.sequence.slice(first, end), qualities: record.qualities.slice(first, end),
    originalLength: record.sequence.length, sourceStart: first + 1, sourceEnd: end,
    trimNote: first || end < record.sequence.length ? `原读段保留 ${first + 1}–${end} nt` : "质量修剪：保留全长" };
}
