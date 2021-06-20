Generates tag density heatmap and histogram for the centered list of features in a headerless regions file.

- If provided regions file is a gene list with the following columns `chrom start end name score strand` set `Gene TSS` as a re-centering criteria.
- If provided regions file is a peak list with the following columns `chrom start end name` set `Peak Center` as a re-centering criteria.

`score` column is always ignored.