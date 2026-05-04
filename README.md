# castanetgraph

`castanetgraph.py` processes Castanet per-probe depth/call output together with precomputed graph-path data, then reports top-supported graph paths and optional depth/consensus artifacts.

## Usage

```bash
python castanetgraph.py \
  --inputgraphdata /path/to/graphdata.parquet \
  --castanetfolder /path/to/castanet_sample_folder \
  --output /path/to/output/sample_prefix
```

## Required arguments

- `--inputgraphdata`: Parquet file containing graph/path summaries.
- `--castanetfolder`: Castanet sample output folder containing `*_depth.csv` plus optional depth/consensus subfolders.
- `--output`: Output prefix used for all generated files.

## Expected Castanet folder layout

At minimum:

- `<castanetfolder>/*_depth.csv`

Optional (used for plots/consensus sequence export if present):

- `<castanetfolder>/Depth_output/*depth_by_pos.csv`
- `<castanetfolder>/consensus_data/<probe_or_block_id>/*_remapped_consensus_sequence.fasta`

## Main outputs

Using `--output /path/prefix`, the script writes:

- `/path/prefix_top_paths.tsv` (created only when no paths remain after filtering)
- `/path/prefix_all_targets_w_top_paths.tsv` (main merged call table)
- `/path/prefix_rMLST_combined.tsv` (when rMLST grouping is detected)
- `/path/prefix_depthplots/` (depth plots when depth files exist)
- `/path/prefix_consensus/` (consensus FASTA files when consensus input exists)

