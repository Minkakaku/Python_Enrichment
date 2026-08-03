# R-Python-scores-enrichment-workflow-template

基因集打分、功能富集、dotplot、箱线图和小提琴图代码素材模板。

## 主要代码

- `sc_score.py`、`score/sc_score.py`：Scanpy 基因集打分代码。
- `score/sc_score_boxplot.R`、`score/sc_score_violin.R`：打分结果绘图。
- `functionenrichment.py`：GO/KEGG 等功能富集。
- `sc_dotplot.py`：单细胞 dotplot。
- `score/report.rmd`：结果报告模板。
- `platform_wrappers/`：原平台启动代码。

## 使用前注意

部分 Python 文件依赖 Streamlit 和内部模块，并含平台临时目录。新项目应提取其中的分析函数，按项目 `results/` 规范输出 PNG/PDF 和支撑表格，不直接运行平台页面代码。

原始副本保存在 `../R-Python-bioinformatics-code-archive/legacy/original/08_scores_enrichment/`。

## 其他人员来源代码

`source_catalog/` 收录 294 份打分、富集和可视化语义副本；`SOURCE_CATALOG.tsv` 保留场景、workflow role、owner 和哈希。正式项目仍应按样本层面统计、PNG/PDF 成对输出及支撑表格规范重新组合代码。
