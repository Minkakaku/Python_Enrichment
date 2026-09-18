# R-Python-scores-enrichment-workflow-template

独立评分与富集分析仓库。每个分析点是一个可直接运行的根目录完整脚本：

- module_score_workflow.py：注释后 h5ad + GMT，计算 Scanpy module score，并按样本×细胞类型比较 CON/EXP。
- gsva_workflow.R：注释后 Seurat RDS + GMT，先汇总样本×细胞类型表达，再运行 GSVA。
- enrichment_workflow.R：差异表达 CSV，分开进行 Up/Down GO/KEGG ORA，并用全部可检验基因运行 preranked GSEA。

## Module score

    python module_score_workflow.py --input annotated.h5ad --gmt signatures.gmt --output results/module_score --sample-col sample --group-col group --celltype-col cell_type

## GSVA

    Rscript gsva_workflow.R --input annotated.rds --gmt signatures.gmt --output results/gsva --sample-col sample --group-col group --celltype-col cell_type

两种评分的组间统计均使用生物学样本汇总，不把单个细胞当作重复。delta_EXP_minus_CON > 0 表示 EXP 高。

## ORA 与 GSEA

    Rscript enrichment_workflow.R --input differential_expression.csv --output results/enrichment --species human --gene-col gene --logfc-col avg_log2FC --pvalue-col p_val --padj-col p_val_adj

ORA 按 adjusted P 和 logFC 阈值拆分 Up/Down，使用全部可映射基因作为背景。GSEA 对全部可检验基因排序，不先按 DEG 显著性过滤。正 NES 表示 EXP 端富集，负 NES 表示 CON 端富集。

评分和富集都是转录层面的相对证据，不直接等于蛋白活性或代谢通量。旧 Streamlit、固定路径和细胞级伪重复脚本原样保存在 legacy/original/。

方法依据：GSVA（PMID: 23323831）、clusterProfiler（PMID: 22455463、34557778）及 Seurat/Scanpy 官方 module-score 文档。