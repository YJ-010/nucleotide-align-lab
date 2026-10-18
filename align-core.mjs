const VALID_BASES = /^[ACGTRYSWKMBDHVN]+$/;
const COMPLEMENT = {
  A: "T", C: "G", G: "C", T: "A", R: "Y", Y: "R", S: "S",
  W: "W", K: "M", M: "K", B: "V", D: "H", H: "D", V: "B", N: "N"
};
const MASK = { A: 1, C: 2, G: 4, T: 8, R: 5, Y: 10, S: 6, W: 9,
  K: 12, M: 3, B: 14, D: 13, H: 11, V: 7, N: 15 };

function normalizeSequence(raw, name) {
  // Coordinate numbers and spacing are common in copied sequence views.
  const sequence = raw.replace(/[\s\d-]/g, "").toUpperCase().replace(/U/g, "T");
  if (!sequence) throw new Error(`“${name}”没有可用的核苷酸序列。`);
  if (!VALID_BASES.test(sequence)) {
    const bad = sequence.match(/[^ACGTRYSWKMBDHVN]/)?.[0];
    throw new Error(`“${name}”包含无效字符“${bad}”；仅支持 DNA 与 IUPAC 模糊碱基。`);
  }
  return sequence;
}

function uniqueNames(records) {
  const count = new Map();
  return records.map((record, i) => {
    const original = (record.name || `序列 ${i + 1}`).trim().slice(0, 90);
    const used = (count.get(original) || 0) + 1;
    count.set(original, used);
    return { name: used === 1 ? original : `${original} (${used})`,
      sequence: normalizeSequence(record.sequence, original) };
  });
}

function parseFasta(input) {
  const records = [];
  let current = null;
  for (const line of input.split(/\r?\n/)) {
    if (line.trim().startsWith(">")) {
      if (current) records.push(current);
      current = { name: line.trim().slice(1).trim() || `序列 ${records.length + 1}`, sequence: "" };
    } else if (current) {
      current.sequence += line;
    } else if (line.trim()) {
      throw new Error("FASTA 第一行应以 > 开头。");
    }
  }
  if (current) records.push(current);
  return records;
}

function parseGenBank(input) {
  const entries = input.split(/(?=^LOCUS\s+)/gm).filter(part => /^LOCUS\s+/m.test(part));
  return entries.map((entry, i) => {
    const locus = entry.match(/^LOCUS\s+(\S+)/m)?.[1];
    const accession = entry.match(/^ACCESSION\s+(\S+)/m)?.[1];
    const originStart = entry.search(/^ORIGIN[^\n]*\r?\n/m);
    const origin = originStart < 0 ? null :
      entry.slice(originStart).replace(/^ORIGIN[^\n]*\r?\n/, "").split(/^\/\/\s*$/m)[0];
    if (!origin) throw new Error(`GenBank 记录 ${locus || i + 1} 缺少 ORIGIN 序列段。`);
    return { name: locus || accession || `GenBank ${i + 1}`, sequence: origin };
  });
}

function parseClustal(input) {
  const map = new Map();
  for (const line of input.split(/\r?\n/).slice(1)) {
    if (!line.trim()) continue;
    if (/^\s/.test(line)) continue; // consensus or annotation row
    const parts = line.trim().split(/\s+/);
    if (parts.length < 2 || !/^[ACGTRYSWKMBDHVNU.\-]+$/i.test(parts[1])) continue;
    map.set(parts[0], (map.get(parts[0]) || "") + parts[1].replace(/\./g, "-"));
  }
  if (!map.size) throw new Error("CLUSTAL 文件中没有找到可读的序列行。");
  return [...map].map(([name, sequence]) => ({ name, sequence }));
}

export function parseSequences(input, fallbackName = "序列 1") {
  const text = input.replace(/^\uFEFF/, "").trim();
  if (!text) return [];
  let records;
  if (/^>/m.test(text) && /^\s*>/.test(text)) records = parseFasta(text);
  else if (/^LOCUS\s+/m.test(text) && /^\s*LOCUS\s+/.test(text)) records = parseGenBank(text);
  else if (/^(CLUSTAL|MUSCLE)\b/i.test(text)) records = parseClustal(text);
  else records = [{ name: fallbackName, sequence: text }];
  return uniqueNames(records);
}

export function reverseComplement(sequence) {
  return [...sequence].reverse().map(base => COMPLEMENT[base] || "N").join("");
}

function pairwise(reference, query) {
  const n = reference.length, m = query.length;
  const stride = m + 1;
  const trace = new Uint8Array((n + 1) * stride);
  let previous = new Int32Array(stride);
  let current = new Int32Array(stride);
  for (let j = 1; j <= m; j++) { previous[j] = -2 * j; trace[j] = 2; }
  for (let i = 1; i <= n; i++) {
    current[0] = -2 * i;
    trace[i * stride] = 1;
    for (let j = 1; j <= m; j++) {
      const a = reference[i - 1], b = query[j - 1];
      const pairScore = a === b ? 2 : (MASK[a] & MASK[b]) ? 0 : -1;
      const diagonal = previous[j - 1] + pairScore;
      const up = previous[j] - 2;
      const left = current[j - 1] - 2;
      let direction = 0, best = diagonal;
      if (up > best) { best = up; direction = 1; }
      if (left > best) { best = left; direction = 2; }
      current[j] = best;
      trace[i * stride + j] = direction;
    }
    [previous, current] = [current, previous];
  }
  const score = previous[m], ref = [], target = [];
  let i = n, j = m;
  while (i || j) {
    const direction = trace[i * stride + j];
    if (direction === 0) { ref.push(reference[--i]); target.push(query[--j]); }
    else if (direction === 1) { ref.push(reference[--i]); target.push("-"); }
    else { ref.push("-"); target.push(query[--j]); }
  }
  return { reference: ref.reverse().join(""), query: target.reverse().join(""), score };
}

function mergeAlignment(rows, incoming) {
  const oldRef = rows[0], newRef = incoming.reference, newQuery = incoming.query;
  const builders = rows.map(() => []);
  const appended = [];
  let i = 0, j = 0;
  while (i < oldRef.length || j < newRef.length) {
    const oldGap = i < oldRef.length && oldRef[i] === "-";
    const newGap = j < newRef.length && newRef[j] === "-";
    if (oldGap && newGap) {
      rows.forEach((row, k) => builders[k].push(row[i]));
      appended.push(newQuery[j]); i++; j++;
    } else if (oldGap || j >= newRef.length) {
      rows.forEach((row, k) => builders[k].push(row[i]));
      appended.push("-"); i++;
    } else if (newGap || i >= oldRef.length) {
      rows.forEach((_, k) => builders[k].push("-"));
      appended.push(newQuery[j]); j++;
    } else {
      if (oldRef[i] !== newRef[j]) throw new Error("序列比对合并时参照序列不一致。");
      rows.forEach((row, k) => builders[k].push(row[i]));
      appended.push(newQuery[j]); i++; j++;
    }
  }
  return [...builders.map(parts => parts.join("")), appended.join("")];
}

function summarize(rows, records, orientations) {
  const length = rows[0].length;
  const consensus = [];
  const allMatch = [];
  for (let i = 0; i < length; i++) {
    const column = rows.map(row => row[i]);
    const exact = column.every(base => base === column[0] && base !== "-");
    consensus.push(exact ? column[0] : "·");
    allMatch.push(exact && /[ACGT]/.test(column[0]) ? "*" : " ");
  }
  const stats = rows.slice(1).map((row, index) => {
    let matches = 0, substitutions = 0, gaps = 0, compared = 0, ambiguous = 0;
    const variants = [];
    let refPosition = 0, queryPosition = 0;
    for (let i = 0; i < length; i++) {
      const a = rows[0][i], b = row[i];
      if (a !== "-") refPosition++;
      if (b !== "-") queryPosition++;
      if (a === "-" || b === "-") {
        gaps++;
        if (variants.length < 80) variants.push({ position: a === "-" ? `${refPosition} 后` : String(refPosition),
          reference: a, alternative: b, kind: "插入/缺失" });
      } else {
        compared++;
        if (a === b && /[ACGT]/.test(a)) matches++;
        else if (a === b || (MASK[a] & MASK[b])) ambiguous++;
        else {
          substitutions++;
          if (variants.length < 80) variants.push({ position: String(refPosition),
            reference: a, alternative: b, kind: "替换" });
        }
      }
    }
    return {
      name: records[index + 1].name, orientation: orientations[index + 1],
      matches, substitutions, gaps, ambiguous, compared,
      identity: compared ? matches / compared * 100 : 0, variants
    };
  });
  return { length, consensus: consensus.join(""), allMatch: allMatch.join(""), stats };
}

export async function alignSequences(records, autoOrient = true, progress = () => {}) {
  if (records.length < 2 || records.length > 8) throw new Error("请提供 2–8 条序列。");
  if (records.some(record => record.sequence.length > 3000))
    throw new Error("单条序列上限为 3,000 nt。长序列建议使用专业桌面软件。");
  const work = records[0].sequence.length * records.slice(1).reduce((sum, r) => sum + r.sequence.length, 0);
  if (work > 18_000_000) throw new Error("本次序列组合过大，请缩短序列或分批比对（计算量上限约 1,800 万格）。");
  let rows = [records[0].sequence];
  const orientations = ["正向"];
  for (let i = 1; i < records.length; i++) {
    progress(i, records.length - 1);
    await new Promise(resolve => setTimeout(resolve, 0)); // let progress paint on small screens
    const forward = pairwise(records[0].sequence, records[i].sequence);
    let choice = forward, orientation = "正向";
    if (autoOrient) {
      const reverse = pairwise(records[0].sequence, reverseComplement(records[i].sequence));
      if (reverse.score > forward.score) { choice = reverse; orientation = "反向互补"; }
    }
    rows = mergeAlignment(rows, choice);
    orientations.push(orientation);
  }
  return { records, rows, orientations, ...summarize(rows, records, orientations) };
}

function safeFastaName(name) { return name.replace(/[\r\n>]/g, " ").trim(); }
function wrap(sequence, width = 80) {
  return sequence.match(new RegExp(`.{1,${width}}`, "g"))?.join("\n") || "";
}
export function toAlignedFasta(result) {
  return result.rows.map((row, index) =>
    `>${safeFastaName(result.records[index].name)}${result.orientations[index] === "反向互补" ? " [reverse-complement]" : ""}\n${wrap(row)}`).join("\n") + "\n";
}
export function toClustal(result) {
  const ids = result.records.map((_, i) => `Seq${i + 1}`.padEnd(9));
  const lines = ["CLUSTAL nucleotide alignment", ""];
  for (let at = 0; at < result.length; at += 60) {
    result.rows.forEach((row, i) => lines.push(`${ids[i]}  ${row.slice(at, at + 60)}`));
    lines.push(`${"".padEnd(11)}${result.allMatch.slice(at, at + 60)}`, "");
  }
  return lines.join("\n") + "\n";
}
export function toText(result) {
  const lines = [
    "核苷酸序列比对  |  全局比对（第一条序列为参照）",
    "得分：匹配 +2，模糊兼容 0，不匹配 -1，缺口 -2",
    "一致率 = 完全相同的 A/C/G/T 位点 ÷ 双方均有碱基的位点",
    ...result.records.map((record, i) => `Seq${i + 1} = ${record.name} (${record.sequence.length} nt, ${result.orientations[i]})`),
    "", ...result.stats.map((stat, i) =>
      `Seq${i + 2} / Seq1: 一致率 ${stat.identity.toFixed(1)}% · 完全相同 ${stat.matches} · 替换 ${stat.substitutions} · 缺口列 ${stat.gaps} · 模糊位点 ${stat.ambiguous}`), ""
  ];
  const coords = result.rows.map(() => 0);
  for (let at = 0; at < result.length; at += 60) {
    const end = Math.min(result.length, at + 60);
    lines.push(`列 ${at + 1}–${end}`);
    result.rows.forEach((row, i) => {
      const start = coords[i] + 1;
      const fragment = row.slice(at, end);
      coords[i] += fragment.replace(/-/g, "").length;
      lines.push(`Seq${i + 1}`.padEnd(7) + String(start).padStart(5) + "  " + fragment + "  " + String(coords[i]));
    });
    lines.push("共识".padEnd(7) + "       " + result.consensus.slice(at, end), "");
  }
  return lines.join("\n") + "\n";
}
