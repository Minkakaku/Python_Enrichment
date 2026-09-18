#!/usr/bin/env python
"""Scanpy module-score 独立入口：h5ad + GMT → 单细胞评分 + 样本级比较。"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

GROUP_COLORS = {"CON": "#72b1b9", "EXP": "#f6ae00"}
FEATURE_COLORS = [
    "#FFEFD5", "#E6E6FA", "#87CEFA", "#6495ED",
    "#4169E1", "#0000CD", "#000080",
]


# 功能：读取标准 GMT，每行第一列为集合名、第二列为描述、后续为基因。
def read_gmt(path: Path) -> dict[str, list[str]]:
    gene_sets: dict[str, list[str]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                gene_sets[fields[0]] = list(dict.fromkeys(fields[2:]))
    return gene_sets


# 功能：从样本名稳健推断 CON/EXP；已有干净分组列时不使用本函数。
def infer_group(sample: str) -> str | None:
    text = str(sample).upper()
    if re.search(r"(^|[_\-.])CON(TROL)?([_\-.]|$)", text):
        return "CON"
    if re.search(r"(^|[_\-.])EXP(ERIMENT)?([_\-.]|$)", text):
        return "EXP"
    return None


# 功能：把图同时保存为 300 dpi PNG 和矢量 PDF。
def save_figure(figure: plt.Figure, stem: Path, width: float, height: float) -> None:
    figure.set_size_inches(width, height)
    figure.tight_layout()
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)


# 功能：用 Scanpy 官方 score_genes 为每个基因集计算匹配背景的 module score。
def calculate_scores(
    adata: ad.AnnData,
    gene_sets: dict[str, list[str]],
    ctrl_size: int,
    n_bins: int,
    seed: int,
) -> tuple[list[str], pd.DataFrame]:
    records = []
    score_columns = []
    for name, genes in gene_sets.items():
        matched = [gene for gene in genes if gene in adata.var_names]
        records.append(
            {"gene_set": name, "input_genes": len(genes), "matched_genes": len(matched)}
        )
        if len(matched) < 3:
            continue
        column = f"score_{name}"
        sc.tl.score_genes(
            adata,
            gene_list=matched,
            score_name=column,
            ctrl_size=ctrl_size,
            n_bins=n_bins,
            random_state=seed,
            use_raw=False,
        )
        score_columns.append(column)
    return score_columns, pd.DataFrame(records)


# 功能：按样本和细胞类型汇总均值，避免把单个细胞当作生物学重复。
def summarize_by_sample(
    adata: ad.AnnData,
    score_columns: list[str],
    sample_col: str,
    group_col: str,
    celltype_col: str,
) -> pd.DataFrame:
    columns = [sample_col, group_col, celltype_col, *score_columns]
    table = adata.obs[columns].copy()
    return (
        table.groupby([sample_col, group_col, celltype_col], observed=True)[score_columns]
        .mean()
        .reset_index()
    )


# 功能：在每个细胞类型内以样本均值进行 CON/EXP Wilcoxon 秩和检验。
def compare_groups(
    summary: pd.DataFrame,
    score_columns: list[str],
    group_col: str,
    celltype_col: str,
) -> pd.DataFrame:
    records = []
    for celltype, subset in summary.groupby(celltype_col, observed=True):
        for score in score_columns:
            con = subset.loc[subset[group_col] == "CON", score].dropna()
            exp = subset.loc[subset[group_col] == "EXP", score].dropna()
            pvalue = mannwhitneyu(exp, con, alternative="two-sided").pvalue if len(con) >= 2 and len(exp) >= 2 else np.nan
            records.append(
                {
                    "cell_type": celltype,
                    "score": score,
                    "n_CON": len(con),
                    "n_EXP": len(exp),
                    "mean_CON": con.mean(),
                    "mean_EXP": exp.mean(),
                    "delta_EXP_minus_CON": exp.mean() - con.mean(),
                    "pvalue": pvalue,
                }
            )
    result = pd.DataFrame(records)
    if not result.empty:
        valid = result["pvalue"].notna()
        result["p_adjust_BH"] = np.nan
        result.loc[valid, "p_adjust_BH"] = multipletests(
            result.loc[valid, "pvalue"], method="fdr_bh"
        )[1]
    return result


# 功能：在固定 UMAP 坐标上绘制每个 module score。
def plot_umap_scores(adata: ad.AnnData, score_columns: list[str], output: Path) -> None:
    if "X_umap" not in adata.obsm:
        return
    cmap = LinearSegmentedColormap.from_list("hongf_feature", FEATURE_COLORS)
    coordinates = adata.obsm["X_umap"]
    for score in score_columns:
        values = adata.obs[score].to_numpy()
        lower, upper = np.nanquantile(values, [0.05, 0.95])
        figure, axis = plt.subplots()
        points = axis.scatter(
            coordinates[:, 0], coordinates[:, 1], c=values, s=2,
            cmap=cmap, vmin=lower, vmax=upper, rasterized=True,
        )
        axis.set_title(score.removeprefix("score_"))
        axis.set_aspect("equal")
        axis.set_axis_off()
        figure.colorbar(points, ax=axis, fraction=0.04)
        save_figure(figure, output / f"UMAP_{score}", 7, 6)


# 功能：按细胞类型展示每个样本的 CON/EXP score 均值。
def plot_sample_scores(
    summary: pd.DataFrame,
    score_columns: list[str],
    group_col: str,
    celltype_col: str,
    output: Path,
) -> None:
    for score in score_columns:
        figure, axis = plt.subplots()
        sns.boxplot(
            data=summary, x=celltype_col, y=score, hue=group_col,
            hue_order=["CON", "EXP"], palette=GROUP_COLORS,
            showfliers=False, ax=axis,
        )
        sns.stripplot(
            data=summary, x=celltype_col, y=score, hue=group_col,
            hue_order=["CON", "EXP"], dodge=True, color="black",
            size=3, alpha=0.75, legend=False, ax=axis,
        )
        axis.tick_params(axis="x", rotation=75)
        axis.set_xlabel("")
        axis.set_ylabel("Sample-level mean score")
        axis.set_title(score.removeprefix("score_"))
        sns.despine(ax=axis)
        save_figure(figure, output / f"sample_level_{score}", 11, 7)


# 功能：定义输入 h5ad、GMT、metadata 列和评分参数。
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Calculate Scanpy module scores from h5ad and GMT.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--gmt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sample-col", default="sample")
    parser.add_argument("--group-col", default="group")
    parser.add_argument("--celltype-col", default="cell_type")
    parser.add_argument("--layer", default=None)
    parser.add_argument("--use-raw", action="store_true")
    parser.add_argument("--ctrl-size", type=int, default=50)
    parser.add_argument("--n-bins", type=int, default=25)
    parser.add_argument("--seed", type=int, default=1001)
    return parser


# 功能：串联评分、样本级汇总、组间比较、绘图和结果说明。
def main() -> None:
    args = build_parser().parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.input)
    if args.use_raw:
        expression = adata.raw.to_adata()
        expression.obs = adata.obs.copy()
        # raw 只保存表达矩阵和基因信息；显式带回 UMAP 等细胞坐标供结果绘图。
        for key, value in adata.obsm.items():
            expression.obsm[key] = value.copy()
        adata = expression
    elif args.layer:
        adata.X = adata.layers[args.layer].copy()
    if args.group_col not in adata.obs:
        adata.obs[args.group_col] = adata.obs[args.sample_col].map(infer_group)
    adata.obs[args.group_col] = pd.Categorical(adata.obs[args.group_col], ["CON", "EXP"])

    score_columns, mapping = calculate_scores(
        adata, read_gmt(args.gmt), args.ctrl_size, args.n_bins, args.seed
    )
    mapping.to_csv(args.output / "gene_set_mapping_summary.csv", index=False)
    adata.write_h5ad(args.output / "module_score_annotated.h5ad")
    summary = summarize_by_sample(
        adata, score_columns, args.sample_col, args.group_col, args.celltype_col
    )
    summary.to_csv(args.output / "sample_level_score_summary.csv", index=False)
    compare_groups(summary, score_columns, args.group_col, args.celltype_col).to_csv(
        args.output / "sample_level_CON_EXP_statistics.csv", index=False
    )
    plot_umap_scores(adata, score_columns, args.output)
    plot_sample_scores(summary, score_columns, args.group_col, args.celltype_col, args.output)

    (args.output / "README_results_CN.md").write_text(
        "\n".join(
            [
                "# Module score 结果说明", "",
                f"- 输入 h5ad：{args.input}", f"- GMT：{args.gmt}",
                f"- 样本列/分组列/细胞类型列：{args.sample_col}/{args.group_col}/{args.celltype_col}",
                "- 单细胞分数由 scanpy.tl.score_genes 计算；组间统计使用样本×细胞类型均值。",
                "- delta_EXP_minus_CON > 0 表示 EXP 较高。",
                "- module score 是相对表达评分，不直接等于通路活性、蛋白活性或代谢通量。",
            ]
        ) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
