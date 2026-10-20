# Random clustering baseline for the omni-pipelines-plan benchmark.
# Assigns every cell to a random cluster, so a real method should score
# well above it. Writes the same outputs as the other methods modules.

import argparse
import json
import os

import numpy as np
import pandas as pd
import scanpy as sc

parser = argparse.ArgumentParser(description="Random clustering baseline")
parser.add_argument("--output_dir", "-o", type=str, required=True, help="Output directory")
parser.add_argument("--name", "-n", type=str, required=True, help="Dataset name")
parser.add_argument("--data.h5ad", dest="data_h5ad", type=str, required=True, help="Input h5ad file")
parser.add_argument(
    "--data.clusters_truth_num",
    dest="clusters_truth_num_path",
    type=str,
    required=True,
    help="File with the true number of clusters",
)
parser.add_argument("--seed", type=int, default=42, help="Seed for the random number generator")
parser.add_argument(
    "--d_cluster",
    type=int,
    default=0,
    help="Number of clusters relative to the true number of clusters",
)
parser.add_argument("--n_comp", type=int, default=30, help="Number of PCs for the structure metrics")
args, _ = parser.parse_known_args()

# data
adata = sc.read_h5ad(args.data_h5ad)
adata.var_names_make_unique()
with open(args.clusters_truth_num_path) as f:
    n_cluster = int(f.read()) + args.d_cluster

# random clusters: each cell gets a cluster from 1..n_cluster
rng = np.random.default_rng(args.seed)
labels = rng.integers(1, n_cluster + 1, size=adata.n_obs)

clusters = pd.DataFrame(
    {"cell_id": adata.obs_names, "louvain": labels, "leiden": labels}
)
clusters.to_csv(
    os.path.join(args.output_dir, f"{args.name}.clusters.tsv"), sep="\t", index=False
)

# PCA of the data, only needed by the metrics module for the structure scores
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
n_comp = min(args.n_comp, adata.n_obs - 1, adata.n_vars - 1)
sc.pp.pca(adata, n_comps=n_comp, random_state=args.seed)
pca_df = pd.DataFrame(
    adata.obsm["X_pca"],
    index=adata.obs_names,
    columns=[f"PC{i + 1}" for i in range(adata.obsm["X_pca"].shape[1])],
)
pca_df.to_csv(
    os.path.join(args.output_dir, f"{args.name}.pca.tsv"),
    sep="\t",
    index=True,
    index_label="cell_id",
)

# timings and clustering info: same structure as the other methods, nothing to report
timings = {
    "gpu_load": None,
    "find_mit_gene": None,
    "filter": None,
    "normalization": None,
    "hvg": None,
    "scaling": None,
    "pca": None,
    "t_sne": None,
    "umap": None,
    "louvain": None,
    "leiden": None,
}
clustering_info = {
    "resolutions": {"louvain": None, "leiden": None},
    "num_runs": {"louvain": None, "leiden": None},
}
with open(os.path.join(args.output_dir, f"{args.name}.timings.json"), "w") as f:
    json.dump(timings, f, indent=2)
with open(os.path.join(args.output_dir, f"{args.name}.clustering_info.json"), "w") as f:
    json.dump(clustering_info, f, indent=2)
