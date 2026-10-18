# SeqLab · 核苷酸序列比对

公开网站：https://yj-010.github.io/nucleotide-align-lab/

免费、免登录的中文 DNA 片段比对工具，计算在浏览器本地完成，序列不上传服务器。采用独立 GitHub Pages 网站。

## 使用

1. 粘贴序列或上传文件；第一条序列作为参照。
2. 点击“开始比对”，查看逐位比对、坐标、一致率与差异概览。
3. 下载对齐 FASTA、CLUSTAL 或纯文本报告。

输入支持 FASTA、GenBank/ApE 文本（.ape/.gb）、CLUSTAL 和纯序列。支持多文件上传、拖放和反向互补识别。结果采用熟悉的分段序列显示方式；不读写 DNAMAN 专有工程文件。

## 范围与算法

- 2–8 条 DNA 序列，单条最多 3,000 nt；动态规划计算量上限为 1,800 万格。
- 两序列采用 Needleman–Wunsch 全局比对及线性缺口罚分：匹配 +2、模糊兼容 0、不匹配 -1、缺口 -2。
- 多序列以第一条为参照合并两两比对；不是完整的渐进式多序列比对流程。
- 适合日常基因片段检查；复杂长片段、结构变异及发表级分析应使用合适的专业方法复核。

## 发布

GitHub Pages：Settings → Pages → Deploy from a branch → main → /(root)。

根目录包含 index.html、styles.css、app.mjs、align-core.mjs 和 .nojekyll，无需构建或外部脚本 CDN。修改 main 分支后自动发布。
