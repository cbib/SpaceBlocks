"""
generate_annotation_template.py
===============================
Write an empty cluster-annotation template TSV for the user to fill in.

Rows    : cluster ids 0 .. (max_clusters - 1)
Columns : "{sample_id}_{resolution}" for every sample × resolution in the scan.

The user enters cell-type labels, then points config["cluster_annotations"] at this
file to unlock annotate_cells / pseudobulk / neighbourhood. Moved out of a `run:`
directive into a script per the Snakemake linter recommendation.
"""

import pandas as pd

params = snakemake.params
columns = [f"{s}_{r}" for s in params.sample_ids for r in params.resolutions]
df = pd.DataFrame(index=range(params.max_clusters), columns=columns, data="")
df.index.name = "cluster"
df.to_csv(str(snakemake.output.template), sep="\t")
