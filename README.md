# SeqLab · 序列比对与克隆设计

公开网站：https://yj-010.github.io/nucleotide-align-lab/

免费、免登录的中文 DNA 片段比对工具，计算在浏览器本地完成，序列不上传服务器。独立部署于 GitHub Pages。

## 测序公司返回的文件怎么用

1. 粘贴或上传目标参照序列（FASTA、GenBank/ApE、SEQ/TXT 等）。
2. 拖入测序公司返回的 `.ab1`、`.seq`、FASTA 或少量 Phred+33 FASTQ 文件。压缩包须先解压，同一样本不必重复导入 AB1 和 SEQ。
3. 确认第一条为目标参照；导入顺序不对时，可在记录卡片中点击“设为参照”。
4. 使用“测序片段定位”，保留自动识别反向互补。需要时勾选两端质量修剪，再点击“开始比对”。
5. 查看参考覆盖区间、方向、一致率、替换和缺口。AB1 可点击“查看峰图”。可导出对齐 FASTA、CLUSTAL 或纯文本报告。

“试试测序片段示例”提供明确标记的模拟数据，用于演示反向读段、质量修剪和单个替换。

## 文件和质量

- AB1/ABI：按 ABIF 签名识别，读取已判读碱基；优先 PBAS2，缺失时使用 PBAS1。质量和峰位置只使用同一套 PCON/PLOC 标签。按 FWO_1 顺序读取 DATA9–12 作为四通道峰图。
- SEQ/TXT、FASTA、GenBank/ApE 文本和 CLUSTAL：自动识别常用序列格式，保留 IUPAC 模糊碱基。
- FASTQ：仅支持明确按 Phred+33 编码的少量读段；不自动区分旧 Phred+64。支持换行的序列/质量，严格检查长度。
- 质量缺失或不一致时不生成虚假质量值。界面显示 Q≥20 比例与模糊碱基数。
- 可选修剪：在读段中找到首个和末个连续 10 个明确 A/C/G/T 且 Q≥20 的区间，保留两者之间全部碱基。只去两端，不去中间、不裁剪参照。无质量值时保留全长；没有合格区间时报错。
- 原始序列保留在输入框，报告记录保留的原读段区间。手动改动记录名称或序列会取消该记录的质量/峰图关联，重新上传可恢复。
- 峰图始终使用原文件方向及原始读段坐标，独立于修剪和反向互补设置。并非重新判读、混合峰检测或杂合突变自动判定。

## 比对范围

- 2–8 条 DNA 序列；测序片段模式参照 ≤12,000 nt、读段 ≤3,000 nt；全长模式每条 ≤3,000 nt。
- 动态规划总计算量上限 1,800 万格；较长参照或多个读段应分批。
- 全局 Needleman–Wunsch，线性缺口罚分：明确匹配 +2、模糊兼容 0、不匹配 -1、缺口 -2。
- 测序片段定位采用 fitting alignment：整条读段参与比对，参照两端不罚分。输出仍保留完整参照；两端未覆盖区间不计作缺失，内部缺口计入差异。多个样本共有的空缺列不产生虚假的差异。
- 多序列以第一条为中心合并，不是完整的渐进式多序列比对流程，不进行双向读段拼接或基因鉴定。重复区可能存在多个同分位置，页面只显示其中一个最优位置。
- 一致率 = 明确相同的 A/C/G/T 位点数 ÷ 双方均有碱基的位点数。必须结合覆盖区间、质量和原始峰图解读。
- 不支持大批量高通量 FASTQ、图片/PDF 报告、SCF 或 DNAMAN 专有工程文件。复杂长片段、结构变异及发表级分析应采用合适的专业方法复核。

## 开发与发布

根目录为静态网站：index.html、styles.css、app.mjs、align-core.mjs、sequencing.mjs 和 .nojekyll，无需构建或外部脚本 CDN。

使用 Node.js 22+ 验证核心逻辑：

```sh
node --test tests/*.test.mjs
```

GitHub Pages：Settings → Pages → Deploy from a branch → main → /(root)。修改 main 分支后自动发布。

2026-10-19 验证：17 项核心/格式测试通过；另用 Biopython 仓库中 7 份有效 AB1 样例核对，统一碱基大小写后，序列、质量和四通道峰图数据均与 Biopython 解析一致；无效文件被拒绝。使用 DOM 模拟测试上传、参照切换、质量修剪、峰图翻页及文本导出。未进行真实浏览器视觉复核或手机实机测试。

格式依据：[Applied Biosystems ABIF 规范](https://archive.gfjc.fiu.edu/workshops/resources/literature/ABIF_File_Format.pdf)、[NCBI FASTQ 格式说明](https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/)、[Biopython ABI 解析器文档](https://biopython.org/docs/latest/api/Bio.SeqIO.AbiIO.html)。公开验证样例来自 [Biopython Tests/Abi](https://github.com/biopython/biopython/tree/master/Tests/Abi)，不随本网站分发。

## 猕猴桃跨版本克隆工作区（2026-10-19）

**当前发布状态：手动序列比较、克隆引物设计及导出已可用。自动基因编号查询所需的完整参考库尚未发布。** 四个版本为 Hong Yang v2、Hong Yang v3、Red5 和 White。可点击“试用模拟序列”查看交互；该演示不代表真实参考基因。

完整库的 GitHub 上传被自动审批拒绝，原因是大型服务器来源数据拟公开发布。当前部署提交仅含界面、算法、测试、建库脚本和待发布目录信息，不含参考序列分片。准备的公开范围见 [REFERENCE_RELEASE.md](REFERENCE_RELEASE.md)。以下建库说明描述已实现但尚未接入线上数据的查询能力。

入口：`#clone`。输入基因或转录本编号 → 查看四个版本的同源候选 → 比较 CDS 或上游区 → 选模板 → 生成并导出引物。可将模板直接送入测序比对工作区。

参考版本：Hong Yang v2、Hong Yang v3 update、Red5 PS1_1.69.0（GenBank GCA_003024255.1）、White v1.0。Red5 的 `CEY00_Acc…` 原始 locus_tag 保留为别名，主显示使用 `Acc…`；转录本仍可按原始标识符精确检索。参考库仅含公开参考数据，不包含用户的实验数据。

数据在 `data/kiwi-20261019/`，按版本 / FNV-1a 哈希分为 64 个 gzip JSON 文件；初次只下载索引及命中的分片。现代浏览器使用原生 `DecompressionStream` 解压，比较和引物初筛在本地 Web Worker 完成。没有 ChatGPT 或登录依赖，没有外部 CDN 脚本。所选基因对应分片会被静态托管服务器请求；自有测序序列和引物输入不上传。

### 数据与同源证据

- `scripts/build_reference.py config.json OUTPUT_DIR --threads 4` 可从只读 FASTA + GFF3 重新建库；配置字段为 `id,label,genome,gff,source,cds(optional),protein(optional)`，依赖 Python 3 标准库和 DIAMOND。源文件 SHA-256、排除模型数、配套 CDS 一致性校验记录在 manifest。
- CDS 按同版本注释从基因组拼接，统一转录方向；仅处理首段 phase，不重复去除跨外显子密码子的相位。重叠 CDS 模型被排除。基因查询默认选最长 CDS；精确转录本编号保留其实际序列。
- 上游区为注释 CDS 5′ 边界前至多 **2000 bp（2 kb）**，排除 CDS；负链反向互补。不是经过实验验证的 TSS 启动子，可能包括 5′UTR、邻近基因或重复序列。不足长度、无 ATG、内部终止、N 和注释异常明确提示。
- 蛋白搜索使用每个基因最长 CDS 翻译产物；DIAMOND blastp `--sensitive --evalue 1e-10 --max-target-seqs 5 --max-hsps 1`，随后保留 identity ≥30%、query/target coverage ≥40% 的前三个候选。互为最佳及同分最佳单独标注；没有把序列相似性当成已证实直系同源，未加入共线性证据。
- HYv3 更新基因组提取 CDS 与较早的配套 CDS 有部分差异，以 manifest 及基因提示为准。White 3 个重叠 CDS 模型（DTZ79_05g01990、DTZ79_15g11550、DTZ79_19g06890）排除。

### 序列比较与引物

- 核苷酸全局 Needleman–Wunsch：+2 / −1 / −2。全长一致率 = A/C/G/T 完全一致列 / 全部对齐列（含缺口及模糊位点）；另外显示无缺口配对覆盖。无自动反向互补，因为参考序列已统一为转录方向。每对最多 3600 万格、单条最多 12,000 nt；超限不伪造结果。
- 克隆引物固定于选定片段两端（默认 18–30 nt、目标 Tm 60°C），不偷偷缩短片段。可显式去掉 CDS 终止密码子，可输入已确定的 5′ 附加序列；产物长度含 / 不含附加序列分别显示。
- Tm 使用 SantaLucia 最近邻、熵盐校正、von Ahsen Mg/dNTP 钠当量；对称校正；已与独立 Biopython `Tm_NN` / `DNA_NN3` / `saltcorr=5` 对照。默认总寡核苷酸 250 nM、单价盐 50 mM、Mg 1.5 mM、总 dNTP 0.2 mM，均可修改。
- 初筛：目标 Tm ±5°C、两端 ΔTm ≤3°C、GC 25–75%、同聚物 ≤5、3′连续互补 <5 bp、简化发夹连续茎 <6 bp；附加序列纳入连续互补检查。仅是计算候选，**不是 Primer3 引擎输出**；没有完整二级结构热力学和全基因组特异性搜索。无法满足固定端点条件时明确返回无候选。
- 对已选参考序列进行引物局部无缺口匹配初筛，报告坐标、错配、末 5 nt 错配及同分位置；不据此宣称全基因组特异。

方法来源：[KGD](https://www.nature.com/articles/s41438-020-0338-9)、[DIAMOND](https://github.com/bbuchfink/diamond)、[Biopython Tm 文档](https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html)。未采用第三方 primd 包，也未打包 Primer3。

### 本次验证

`node --test tests/*.test.mjs` 包含原有测序功能及新的克隆坐标、长度限制、对齐分母、独立 Tm 对照、互补方向、引物端点和导出相关核心验证。另执行 DOM 流程测试（真实参考查询、四版本比较、模板切换、引物设计及送回测序比对）。当前环境未完成真实浏览器视觉与手机实机检查。

手动输入区仅接受每格一条 FASTA 或纯 DNA：CDS 60–12,000 bp，上游区 60–2000 bp。没有坐标时显示“未核对”，不生成虚构染色体位置；支持单条设计引物、两至四个版本比较，编辑后清除旧结果。
