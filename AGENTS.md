# 评分与富集项目约定

- 本仓库是独立分析点；每个根目录脚本接收 RDS、H5AD、GMT 或差异表。
- 所有 CON/EXP 统计在样本层面完成，不把细胞当生物学重复。
- CON 使用 #72b1b9，EXP 使用 #f6ae00；delta 定义为 EXP-CON。
- Up/Down 分开 ORA；GSEA 使用全部可检验基因排序。
- 图成对保存 PNG/PDF，绘图数据和 README_results_CN.md 放同一输出目录。
- legacy/original/ 和 source_catalog/ 保持原字节。