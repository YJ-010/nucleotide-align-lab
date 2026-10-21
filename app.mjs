import { alignSequences, reverseComplement, toAlignedFasta, toClustal, toText } from "./align-core.mjs?v=20261019-clone4";
import { parseInputText, readSequenceFile, qualitySummary, trimRead } from "./sequencing.mjs?v=20261019-clone4";

const $ = id => document.getElementById(id);
const input = $("sequence-input");
const status = $("status");
const results = $("results");
const resultState = $("result-state");
const alignButton = $("align-button");
let currentResult = null;
let view = "difference";
let importedRecords = [];
let inputRevision = 0;
let importing = false;
let traceRecord = null, traceStart = 0;

function invalidateResults() {
  currentResult = null;
  results.hidden = true;
  $("empty-state").hidden = false;
  resultState.textContent = "等待比对";
  resultState.classList.remove("ready");
}

function currentRecords() {
  return parseInputText(input.value).map(record => {
    const saved = importedRecords.find(old => old.name === record.name && old.sequence === record.sequence);
    return saved ? { ...saved } : record;
  });
}

function setRecords(records) {
  const used = new Set();
  importedRecords = records.map(record => {
    const root = record.name.replace(/[\r\n>]/g, " ").trim().slice(0, 80) || "序列";
    let name = root, suffix = 2;
    while (used.has(name)) name = `${root} (${suffix++})`;
    used.add(name);
    return { ...record, name };
  });
  input.value = importedRecords.map(record => `>${record.name}\n${record.sequence}`).join("\n");
  inputRevision++;
  invalidateResults();
  $("trace-panel").hidden = true;
  traceRecord = null;
  renderReadList(importedRecords);
}

function renderReadList(records) {
  const container = $("read-list");
  container.replaceChildren();
  records.forEach((record, index) => {
    const summary = qualitySummary(record);
    const card = document.createElement("div");
    card.className = "read-card";
    const name = document.createElement("strong");
    name.textContent = `${index === 0 ? "参照" : "读段 " + index} · ${record.name}`;
    const details = document.createElement("span");
    details.className = "read-details";
    details.textContent = `${record.format || "文本"} · ${summary.length.toLocaleString()} nt · 模糊碱基 ${summary.ambiguous} · ${summary.q20 === null ? "无质量值" : "Q≥20：" + summary.q20.toFixed(1) + "%"}`;
    card.append(name, details);
    const actions = document.createElement("div");
    actions.className = "read-actions";
    if (index > 0) {
      const referenceButton = document.createElement("button");
      referenceButton.type = "button";
      referenceButton.className = "text-button";
      referenceButton.textContent = "设为参照";
      referenceButton.disabled = importing || alignButton.disabled;
      referenceButton.addEventListener("click", () => {
        if (importing || alignButton.disabled) return;
        const current = currentRecords();
        const selected = current.splice(index, 1)[0];
        setRecords([selected, ...current]);
        setStatus(`已将“${selected.name}”设为参照。`);
      });
      actions.append(referenceButton);
    }
    if (record.trace) {
      const traceButton = document.createElement("button");
      traceButton.type = "button";
      traceButton.className = "text-button";
      traceButton.textContent = "查看峰图";
      traceButton.addEventListener("click", () => {
        traceRecord = record; traceStart = 0; renderTrace();
        $("trace-panel").scrollIntoView({ behavior: "smooth", block: "start" });
      });
      actions.append(traceButton);
    }
    card.append(actions);
    for (const warning of record.warnings || []) {
      const notice = document.createElement("small");
      notice.className = "read-warning";
      notice.textContent = warning;
      card.append(notice);
    }
    container.append(card);
  });
}

function renderTrace() {
  if (!traceRecord?.trace) return;
  const record = traceRecord, { positions, channels } = record.trace;
  const end = Math.min(traceStart + 60, positions.length);
  const firstSample = Math.max(0, positions[traceStart] - 10);
  const lastSample = Math.min(channels.A.length - 1, positions[end - 1] + 10);
  const width = 1120, top = 58, bottom = 230;
  const x = sample => 20 + 1080 * (sample - firstSample) / Math.max(1, lastSample - firstSample);
  let peak = 1;
  for (const channel of Object.values(channels))
    for (let i = firstSample; i <= lastSample; i++) peak = Math.max(peak, channel[i]);
  const svg = $("trace-svg");
  svg.replaceChildren();
  svg.setAttribute("viewBox", `0 0 ${width} 250`);
  function element(tag, attrs, value) {
    const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
    for (const [key, val] of Object.entries(attrs)) node.setAttribute(key, String(val));
    if (value !== undefined) node.textContent = value;
    return node;
  }
  for (let i = traceStart; i < end; i++) {
    const at = x(positions[i]);
    if (record.qualities && record.qualities[i] < 20)
      svg.append(element("rect", { x: at - 6, y: 8, width: 12, height: 225, fill: "#f9dfd6", opacity: .7 }));
    const label = element("text", { x: at, y: 38, "text-anchor": "middle", "font-size": 13, fill: "#244838" }, record.sequence[i]);
    label.append(element("title", {}, `碱基 ${i + 1}：${record.sequence[i]}；${record.qualities ? "Q=" + record.qualities[i] : "无质量值"}`));
    svg.append(label);
    if (i === traceStart || (i + 1) % 10 === 0)
      svg.append(element("text", { x: at, y: 18, "text-anchor": "middle", "font-size": 10, fill: "#75867a" }, String(i + 1)));
  }
  for (const [base, color] of Object.entries({ A: "#25964c", C: "#347dcc", G: "#353d48", T: "#cb4b59" })) {
    let path = "";
    for (let i = firstSample; i <= lastSample; i++)
      path += `${i === firstSample ? "M" : "L"}${x(i).toFixed(1)},${(bottom - (bottom - top) * Math.max(0, channels[base][i]) / peak).toFixed(1)} `;
    svg.append(element("path", { d: path, fill: "none", stroke: color, "stroke-width": 1.4 }));
  }
  $("trace-panel").hidden = false;
  $("trace-title").textContent = record.name;
  $("trace-note").textContent = "显示 AB1 文件内已有的四通道信号与碱基判读。这里始终为原始读段、原始方向，不受修剪和反向互补设置影响。";
  $("trace-range").textContent = `原读段 ${traceStart + 1}–${end} / ${positions.length} nt`;
  $("trace-previous").disabled = traceStart === 0;
  $("trace-next").disabled = end >= positions.length;
}

function updateModeNote() {
  $("mode-note").textContent = $("alignment-mode").value === "fitting"
    ? "适合 Sanger 短读段：参照中未测到的两端不计作缺失；内部缺口仍会显示。"
    : "适合比较相近长度的完整序列；全长两端的差异也会计入缺口。";
}

function setStatus(message, error = false) {
  status.textContent = message;
  status.classList.toggle("error", error);
}
function showExample() {
  const reference = "ATGACCTTCGAGGTCATCGTACAGGCTGATCGATGCCAATCGTACCTGAGATCTACGTGGATGATCCGTTACGAGTCCAAGTTGACCTAG";
  const variant = reference.slice(0, 28) + "A" + reference.slice(29, 56) + "T" + reference.slice(56, 71) + reference.slice(73);
  const sample3 = reference.slice(0, 41) + "G" + reference.slice(41, 63) + "C" + reference.slice(64);
  setRecords([{ name: "参照序列 | 野生型", sequence: reference },
    { name: "样本 A | 替换与缺失", sequence: variant },
    { name: "样本 B | 反向互补示例", sequence: reverseComplement(sample3) }]);
  $("alignment-mode").value = "global";
  $("trim-quality").checked = false;
  updateModeNote();
  setStatus("示例已载入。点击“开始比对”查看结果。");
  input.focus();
}

function showSequencingExample() {
  let seed = 17;
  const reference = Array.from({ length: 320 }, () => {
    seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
    return "ACGT"[(seed >>> 24) % 4];
  }).join("");
  const fragment = reference.slice(67, 258);
  const variant = fragment.slice(0, 70) + ({ A: "C", C: "G", G: "T", T: "A" }[fragment[70]]) + fragment.slice(71);
  const read = "NNNNN" + reverseComplement(variant) + "NNNNN";
  setRecords([{ name: "目标参照（模拟）", sequence: reference, format: "示例" },
    { name: "反向测序读段（模拟）", sequence: read, format: "模拟 FASTQ", warnings: [],
      qualities: Array.from(read, base => base === "N" ? 5 : 35) }]);
  $("alignment-mode").value = "fitting";
  $("trim-quality").checked = true;
  updateModeNote();
  setStatus("模拟示例已载入：反向读段、低质量两端和 1 处替换。点击“开始比对”。");
}

function metric(label, number, unit = "") {
  const card = document.createElement("div");
  card.className = "metric";
  const caption = document.createElement("span");
  caption.className = "metric-label";
  caption.textContent = label;
  const value = document.createElement("strong");
  value.textContent = number;
  if (unit) {
    const small = document.createElement("small");
    small.textContent = " " + unit;
    value.append(small);
  }
  card.append(caption, value);
  return card;
}

function classifyBase(result, rowIndex, column) {
  const here = result.rows[rowIndex][column];
  const uncovered = index => result.mode === "fitting" && index > 0 &&
    (column < result.coverageBounds[index].first || column > result.coverageBounds[index].last);
  if (uncovered(rowIndex)) return "uncovered";
  const other = rowIndex === 0
    ? result.rows.slice(1).filter((_, index) => !uncovered(index + 1)).map(row => row[column])
    : [result.rows[0][column]];
  if (!other.length) return "uncovered";
  if (here === "-" || other.includes("-")) return "gap";
  if (other.every(base => base === here && /[ACGT]/.test(here))) return "same";
  if (/[RYSWKMBDHVN]/.test(here) || other.some(base => /[RYSWKMBDHVN]/.test(base))) return "uncertain";
  return "substitution";
}

function makeBaseFragment(result, rowIndex, from, to) {
  const fragment = document.createElement("span");
  fragment.className = "sequence-fragment";
  for (let column = from; column < to; column++) {
    const base = document.createElement("span");
    base.className = "base " + classifyBase(result, rowIndex, column);
    if ((column + 1) % 10 === 0 && column !== to - 1) base.classList.add("group-end");
    const letter = result.rows[rowIndex][column];
    base.dataset.base = letter;
    base.textContent = letter;
    if (base.classList.contains("uncovered")) base.title = "参照两端未覆盖，不计作缺失";
    fragment.append(base);
  }
  return fragment;
}

function renderAlignment(result) {
  const container = $("alignment-view");
  container.replaceChildren();
  container.className = "alignment-view " + (view === "base" ? "base-mode" : "difference");
  const counters = result.rows.map(() => 0);
  for (let from = 0; from < result.length; from += 60) {
    const to = Math.min(result.length, from + 60);
    const block = document.createElement("section");
    block.className = "align-block";
    const top = document.createElement("div");
    top.className = "block-top";
    const left = document.createElement("span");
    left.textContent = `列 ${from + 1}—${to}`;
    const right = document.createElement("span");
    right.textContent = `${Math.ceil(to / 60)} / ${Math.ceil(result.length / 60)}`;
    top.append(left, right);
    block.append(top);
    result.rows.forEach((row, index) => {
      const piece = row.slice(from, to);
      const bases = piece.replace(/-/g, "").length;
      const line = document.createElement("div");
      line.className = "align-row";
      const label = document.createElement("span");
      label.className = "row-label";
      label.textContent = `Seq${index + 1}`;
      label.title = result.records[index].name;
      const start = document.createElement("span");
      start.className = "row-coordinate";
      start.textContent = bases ? String(counters[index] + 1) : "·";
      const end = document.createElement("span");
      end.className = "row-coordinate end";
      counters[index] += bases;
      end.textContent = bases ? String(counters[index]) : "·";
      line.append(label, start, makeBaseFragment(result, index, from, to), end);
      block.append(line);
    });
    const consensus = document.createElement("div");
    consensus.className = "align-row consensus-row";
    const label = document.createElement("span");
    label.className = "row-label";
    label.textContent = "共识";
    const space = document.createElement("span");
    space.className = "row-coordinate";
    const letters = document.createElement("span");
    letters.className = "sequence-fragment";
    for (let i = from; i < to; i++) {
      const character = document.createElement("span");
      character.className = "base";
      if ((i + 1) % 10 === 0 && i !== to - 1) character.classList.add("group-end");
      character.textContent = result.consensus[i];
      letters.append(character);
    }
    consensus.append(label, space, letters);
    block.append(consensus);
    container.append(block);
  }
}

function renderVariants(result) {
  const list = $("variant-list");
  list.replaceChildren();
  const total = result.stats.reduce((sum, stat) => sum + stat.substitutions + stat.gaps, 0);
  $("variant-count").textContent = `共 ${total} 个差异列（按样本累计）`;
  if (!total) {
    const empty = document.createElement("p");
    empty.className = "variant-empty";
    empty.textContent = "未发现明确的替换或缺口位点。";
    list.append(empty);
    return;
  }
  let shown = 0;
  result.stats.forEach((stat, i) => {
    for (const variant of stat.variants) {
      if (shown >= 48) break;
      const item = document.createElement("div");
      item.className = "variant-item";
      const position = document.createElement("span");
      position.className = "variant-position";
      position.textContent = `#${variant.position}`;
      const change = document.createElement("strong");
      change.textContent = `${variant.reference} → ${variant.alternative}`;
      const name = document.createElement("span");
      name.textContent = `Seq${i + 2} · ${variant.kind}`;
      item.append(position, change, name);
      list.append(item);
      shown++;
    }
  });
  if (shown < total) {
    const extra = document.createElement("p");
    extra.className = "variant-empty";
    extra.textContent = `此处展示前 ${shown} 个差异；完整序列请查看上方比对或下载文件。`;
    list.append(extra);
  }
}

function renderResults(result) {
  $("empty-state").hidden = true;
  results.hidden = false;
  resultState.textContent = "比对完成";
  resultState.classList.add("ready");
  const average = result.stats.reduce((sum, stat) => sum + stat.identity, 0) / result.stats.length;
  $("summary-cards").replaceChildren(
    metric("序列数量", String(result.rows.length), "条"),
    metric("对齐长度", result.length.toLocaleString(), "列"),
    metric("平均一致率", average.toFixed(1) + "%")
  );
  $("column-count").textContent = `${result.length.toLocaleString()} COLUMNS`;
  const key = $("sequence-key");
  key.replaceChildren();
  result.records.forEach((record, index) => {
    const entry = document.createElement("span");
    const badge = document.createElement("b");
    badge.textContent = String(index + 1);
    entry.append(badge, document.createTextNode(record.name));
    if (index > 0 && result.orientations[index] === "反向互补") {
      const note = document.createElement("em");
      note.textContent = "反向互补";
      entry.append(note);
    }
    if (record.trimNote) entry.append(document.createTextNode(` · ${record.trimNote}`));
    key.append(entry);
  });
  const coverage = $("coverage-summary");
  coverage.replaceChildren();
  for (const [index, stat] of result.stats.entries()) {
    const line = document.createElement("p");
    line.textContent = `Seq${index + 2} · 一致率 ${stat.identity.toFixed(1)}% · ${stat.matches} 个明确匹配 · ${stat.substitutions} 处替换 · ${stat.ambiguous} 个模糊位点` +
      (result.mode === "fitting" ? ` · 参照覆盖 ${stat.refStart ?? "—"}–${stat.refEnd ?? "—"} nt（${stat.coverage.toFixed(1)}%），未覆盖 ${stat.uncovered} nt` : "");
    coverage.append(line);
    if (stat.identity < 80 || stat.compared < 30) {
      const note = document.createElement("p");
      note.className = "read-warning";
      note.textContent = `Seq${index + 2} 的匹配比例较低或重叠较短，请核对参照与测序质量；此结果不能单独判定样本正确。`;
      coverage.append(note);
    }
  }
  renderAlignment(result);
  renderVariants(result);
}

async function runAlignment() {
  if (alignButton.disabled || importing) return;
  const revision = inputRevision;
  try {
    invalidateResults();
    const rawRecords = currentRecords();
    const records = rawRecords.map((record, index) => index > 0 && $("trim-quality").checked ? trimRead(record) : record);
    alignButton.disabled = true;
    $("result-panel").setAttribute("aria-busy", "true");
    setStatus(`已识别 ${records.length} 条序列，准备比对…`);
    const result = await alignSequences(records, $("auto-orient").checked, (current, total) =>
      setStatus(`正在比对第 ${current} / ${total} 条样本…`), { mode: $("alignment-mode").value });
    if (revision !== inputRevision) throw new Error("输入或设置已变化，请点击“开始比对”重新计算。");
    currentResult = result;
    renderResults(result);
    setStatus(`完成：${records.length} 条序列，${result.length.toLocaleString()} 个对齐列。序列未离开浏览器。`);
    if (window.matchMedia("(max-width: 1100px)").matches)
      $("result-panel").scrollIntoView({ behavior: "smooth", block: "start" });
  } catch (error) {
    setStatus(error.message || "比对失败，请检查输入格式。", true);
  } finally {
    alignButton.disabled = false;
    $("result-panel").removeAttribute("aria-busy");
  }
}

async function addFiles(files) {
  if (!files.length || importing || alignButton.disabled) return;
  const revision = inputRevision;
  importing = true;
  $("upload-button").disabled = true;
  try {
    const records = currentRecords();
    let hasSequencingFile = false;
    for (const file of files) {
      const loaded = await readSequenceFile(file);
      if (!loaded.length) throw new Error(`“${file.name}”没有可用的序列。`);
      records.push(...loaded);
      hasSequencingFile ||= loaded.some(record => record.format === "AB1" || record.format === "FASTQ");
      if (records.length > 8) throw new Error("加入这些文件后超过 8 条序列。本批文件尚未导入，请减少文件数量。");
    }
    if (revision !== inputRevision) throw new Error("读取文件期间输入发生变化，本批文件尚未导入，请重新选择文件。");
    importing = false;
    setRecords(records);
    if (hasSequencingFile) { $("alignment-mode").value = "fitting"; updateModeNote(); }
    setStatus(`已载入 ${files.length} 个文件，共 ${records.length} 条序列。${records.length === 1 ? "请再加入目标参照或另一条序列，并确认第一条为参照。" : "请确认第一条为参照，再点击“开始比对”。"}`);
  } catch (error) {
    setStatus(error.message || "读取文件失败。", true);
  } finally {
    importing = false;
    $("upload-button").disabled = false;
  }
}

function exportValue(format) {
  if (format === "clustal") return { body: toClustal(currentResult), extension: "aln" };
  if (format === "text") return { body: toText(currentResult), extension: "txt" };
  return { body: toAlignedFasta(currentResult), extension: "fasta" };
}

$("example-button").addEventListener("click", showExample);
$("empty-example-button").addEventListener("click", showExample);
$("sanger-example-button").addEventListener("click", showSequencingExample);
$("trace-previous").addEventListener("click", () => { traceStart = Math.max(0, traceStart - 60); renderTrace(); });
$("trace-next").addEventListener("click", () => { if (traceRecord && traceStart + 60 < traceRecord.sequence.length) { traceStart += 60; renderTrace(); } });
$("close-trace").addEventListener("click", () => { $("trace-panel").hidden = true; });
for (const id of ["alignment-mode", "trim-quality", "auto-orient"]) {
  $(id).addEventListener("change", () => {
    inputRevision++;
    invalidateResults();
    updateModeNote();
    setStatus("设置已更新，请点击“开始比对”重新计算。");
  });
}
input.addEventListener("input", () => {
  inputRevision++;
  invalidateResults();
  $("trace-panel").hidden = true;
  traceRecord = null;
  try {
    const records = currentRecords();
    importedRecords = records;
    renderReadList(records);
  } catch { $("read-list").replaceChildren(); }
  setStatus("输入已更新，请点击“开始比对”。手动更改的记录不沿用原始质量或峰图。");
});
$("align-button").addEventListener("click", runAlignment);
$("clear-button").addEventListener("click", () => {
  setRecords([]);
  $("file-input").value = "";
  resultState.textContent = "等待序列";
  setStatus("");
  input.focus();
});
$("upload-button").addEventListener("click", () => $("file-input").click());
$("file-input").addEventListener("change", async event => {
  await addFiles(event.target.files);
  event.target.value = "";
});
const drop = $("drop-zone");
["dragenter", "dragover"].forEach(eventName => drop.addEventListener(eventName, event => {
  event.preventDefault();
  drop.classList.add("drag-over");
}));
["dragleave", "drop"].forEach(eventName => drop.addEventListener(eventName, event => {
  event.preventDefault();
  drop.classList.remove("drag-over");
}));
drop.addEventListener("drop", event => addFiles(event.dataTransfer.files));
document.querySelectorAll(".view-button").forEach(button => button.addEventListener("click", () => {
  view = button.dataset.view;
  document.querySelectorAll(".view-button").forEach(item => {
    item.classList.toggle("active", item === button);
    item.setAttribute("aria-pressed", String(item === button));
  });
  if (currentResult) $("alignment-view").className = "alignment-view " + (view === "base" ? "base-mode" : "difference");
}));
$("copy-button").addEventListener("click", async () => {
  if (!currentResult) return;
  try {
    await navigator.clipboard.writeText(toText(currentResult));
    setStatus("完整比对文本已复制。");
  } catch {
    setStatus("复制未成功，请使用“纯文本报告”下载。", true);
  }
});
$("download-button").addEventListener("click", () => {
  if (!currentResult) return;
  const { body, extension } = exportValue($("export-format").value);
  const url = URL.createObjectURL(new Blob([body], { type: "text/plain;charset=utf-8" }));
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = `KiwiGene_sequence-alignment.${extension}`;
  document.body.append(anchor);
  anchor.click();
  anchor.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
  setStatus(`已下载 ${anchor.download}。`);
});

// Receive a selected cloning template without sending sequence data to a server.
window.addEventListener('seqlab-reference', event => {
  const record = event.detail;
  if (!record || typeof record.sequence !== 'string' || !/^[ACGTRYSWKMBDHVN]+$/.test(record.sequence)) return;
  setRecords([{name:record.name || '目标片段', sequence:record.sequence, format:'目标片段参照'}]);
  $('alignment-mode').value = 'fitting';
  $('alignment-mode').dispatchEvent(new Event('change'));
  setStatus('目标片段已设为参照，再上传 AB1 / SEQ 文件即可比对。仅检查读段覆盖区域，未覆盖位置及载体连接处仍需验证。');
});

