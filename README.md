# SeqLab · 核苷酸序列比对

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
