# Random clustering baseline

A baseline method for the `methods` stage of the omni-pipelines-plan benchmark.
Every cell is assigned to a random cluster (as many clusters as the true number
plus `d_cluster`), using a fixed seed so the result is reproducible.

Parameters:

- `seed`: seed of the random number generator (default 42 and 173)
- `d_cluster`: number of clusters relative to the true number (default 0)
- `n_comp`: number of principal components written to `pca.tsv` for the structure metrics (default 30)

Outputs: `{name}.clusters.tsv`, `{name}.pca.tsv`, `{name}.timings.json`, `{name}.clustering_info.json`.
