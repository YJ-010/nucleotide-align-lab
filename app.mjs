import { parseSequences, alignSequences, reverseComplement, toAlignedFasta, toClustal, toText } from "./align-core.mjs";

const $ = id => document.getElementById(id);
const input = $("sequence-input");
const status = $("status");
const results = $("results");
const resultState = $("result-state");
const alignButton = $("align-button");
let currentResult = null;
let view = "difference";

function setStatus(message, error = false) {
  status.textContent = message;
  status.classList.toggle("error", error);
}
function showExample() {
  const reference = "ATGACCTTCGAGGTCATCGTACAGGCTGATCGATGCCAATCGTACCTGAGATCTACGTGGATGATCCGTTACGAGTCCAAGTTGACCTAG";
  const variant = reference.slice(0, 28) + "A" + reference.slice(29, 56) + "T" + reference.slice(56, 71) + reference.slice(73);
  const sample3 = reference.slice(0, 41) + "G" + reference.slice(41, 63) + "C" + reference.slice(64);
  input.value = `>参照序列 | 野生型\n${reference}\n>样本 A | 替换与缺失\n${variant}\n>样本 B | 反向互补示例\n${reverseComplement(sample3)}`;
  setStatus("示例已载入。点击“开始比对”查看结果。");
  input.focus();
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
  const other = rowIndex === 0
    ? result.rows.slice(1).map(row => row[column])
    : [result.rows[0][column]];
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
    key.append(entry);
  });
  renderAlignment(result);
  renderVariants(result);
}

async function runAlignment() {
  if (alignButton.disabled) return;
  try {
    const records = parseSequences(input.value);
    alignButton.disabled = true;
    $("result-panel").setAttribute("aria-busy", "true");
    setStatus(`已识别 ${records.length} 条序列，准备比对…`);
    const result = await alignSequences(records, $("auto-orient").checked, (current, total) =>
      setStatus(`正在比对第 ${current} / ${total} 条样本…`));
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
  if (!files.length) return;
  try {
    const records = input.value.trim() ? parseSequences(input.value) : [];
    for (const file of files) {
      if (file.size > 2_000_000) throw new Error(`文件“${file.name}”超过 2 MB，请截取所需片段。`);
      const body = await file.text();
      if (body.includes("\0")) throw new Error(`“${file.name}”是二进制文件，请先在原软件中导出 FASTA 或 GenBank 文本。`);
      const filename = file.name.replace(/\.[^.]+$/, "") || file.name;
      records.push(...parseSequences(body, filename));
    }
    input.value = records.map(record => `>${record.name}\n${record.sequence}`).join("\n");
    setStatus(`已载入 ${files.length} 个文件，共 ${records.length} 条序列。`);
  } catch (error) {
    setStatus(error.message || "读取文件失败。", true);
  }
}

function exportValue(format) {
  if (format === "clustal") return { body: toClustal(currentResult), extension: "aln" };
  if (format === "text") return { body: toText(currentResult), extension: "txt" };
  return { body: toAlignedFasta(currentResult), extension: "fasta" };
}

$("example-button").addEventListener("click", showExample);
$("empty-example-button").addEventListener("click", showExample);
$("align-button").addEventListener("click", runAlignment);
$("clear-button").addEventListener("click", () => {
  input.value = "";
  $("file-input").value = "";
  currentResult = null;
  results.hidden = true;
  $("empty-state").hidden = false;
  resultState.textContent = "等待序列";
  resultState.classList.remove("ready");
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
  anchor.download = `sequence-alignment.${extension}`;
  document.body.append(anchor);
  anchor.click();
  anchor.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
  setStatus(`已下载 ${anchor.download}。`);
});
