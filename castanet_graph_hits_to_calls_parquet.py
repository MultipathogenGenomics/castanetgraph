import shutil
import statistics
import glob
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import re
import networkx as nx
from itertools import combinations
import os
from collections import Counter
from Bio import SeqIO, SeqRecord,Seq
from rmlst_stats import rmlst_stats

"""
castanet_graph_hits_to_calls.py

Read a precomputed graph-data parquet (produced by generate_graphdata_parquet.py)
and a Castanet output folder for a sample, then map coverage/read-count metrics
from the Castanet probe-level file onto path/block summaries from the
pangenome graphs. The script filters candidate paths by coverage and read
support, removes subset paths when a superset path has better support,
and selects a top path per (component, graph_name) combination. also generate depth plots and consensus sequences


"""

#TODO : Update to evaluate paths with shared blocks, e.g. A-B-C and A-B-D where A-B is shared check coverage of C and D to assign read proportions to the two paths. (problem is that can't recalculate ncov_mindepth2 easily without per-base coverage)
# in practice probably need something more probabilistic than above for more complex cases
def get_args():
    """
    Parse command-line arguments.

    Returns a namespace with: inputgraphdata (glob/pickle), inputcastanet (csv), output (prefix).
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputgraphdata", type=str, required=True,
                        help="Path to the input graph data pickle file.")
    parser.add_argument("--castanetfolder", type=str, required=True,
                        help="Path to the input Castanet depth CSV file.")
    parser.add_argument("--output", type=str, required=True,
                        help="Prefix for output files")
    args = parser.parse_args()

    args.inputcastanet = glob.glob(args.castanetfolder + "/*_depth.csv")[0]
    args.inputdepthfolder = args.castanetfolder + "/Depth_output"
    args.consensus = args.castanetfolder + "/consensus_data"

    return args

def extract_species_set(graph_name):
    if isinstance(graph_name, str) and graph_name.startswith("bac000"):
        specset = graph_name.split("block")[-1].replace("bac","BACT")

    elif isinstance(graph_name, str) and "block" in graph_name:
        specset = graph_name.split("block")[-1]

    elif isinstance(graph_name, list):
        specset = ", ".join(graph_name)

    elif isinstance(graph_name, set):
        specset = ", ".join(list(graph_name))

    else:
        specset = graph_name
    if specset == "":
        ...
    if specset.startswith("graph"):
        specset = specset[5:]
    return specset

def update_to_longblockordernames(row):
    blockid = list(row["block_id_order"])
    graphname = row["graph_name"].replace("_","").replace("graph","")
    graphname = graphname.replace("BACT00","bac00").replace("16S","16s").replace("23S","23s")
    outid = [f"graph{graphname}block{b}" for b in blockid]
    return outid
def update_to_longblocksetnames(row):
    blockid = list(row["block_id_set"])
    graphname = row["graph_name"].replace("_","").replace("graph","")
    graphname = graphname.replace("BACT00","bac00").replace("16S","16s").replace("23S","23s")
    outid = [f"graph{graphname}block{b}" for b in blockid]
    return outid

def update_species_set(row):
    graphtochange = ["fungimito_graph","fungiamr_graph"]
    if row["graph_name"] in graphtochange:
        path = list(row["pathnames"])
        pathls = list(set([x.split("_")[0] for x in path]))
        pathls = ", ".join(pathls)
        return pathls
    else:
        path = list(row["pathnames"])
        pathls = list(set(["~".join(x.split("_")[0].split("~")[1:]) for x in path]))
        pathls = ", ".join(pathls)
        return pathls


def process_strict_unique_bact_groups(df):
    # 1. Isolate target rows
    mask = df['graph_name'].str.startswith('bac000', na=False)
    working_df = df[mask].copy()
    # Always return a tuple (df, group_unique_map). If there are no working rows,
    # return the original dataframe and an empty metadata map.
    if working_df.empty:
        return df, {}

    # 2. Build graph (Edges = shared species)
    G = nx.Graph()
    species_to_rows = {}
    for idx, row in working_df.iterrows():
        for s in row['species_set']:
            species_to_rows.setdefault(s, []).append(idx)
    for rows in species_to_rows.values():
        for u, v in combinations(rows, 2):
            G.add_edge(u, v)

    # 3. Find Cliques + Locus Constraint
    raw_groups = []
    for clique in nx.find_cliques(G):
        sub_df = working_df.loc[list(clique)].drop_duplicates(subset='graph_name')
        common_s = set.intersection(*sub_df['species_set'].map(set))
        if len(common_s) > 0:
            raw_groups.append({'indices': set(sub_df.index), 'common': list(common_s)})

    # 4. Redundancy Filter
    raw_groups.sort(key=lambda x: len(x['indices']), reverse=True)
    final_valid_groups = []
    for i, g_a in enumerate(raw_groups):
        if not any(i != j and g_a['indices'].issubset(g_b['indices']) for j, g_b in enumerate(raw_groups)):
            final_valid_groups.append(g_a)

    # --- UNIQUE LOCUS CALCULATION ---
    all_indices = [idx for g in final_valid_groups for idx in g['indices']]
    index_counts = Counter(all_indices)
    for g in final_valid_groups:
        g['n_unique_loci'] = sum(1 for idx in g['indices'] if index_counts[idx] == 1)

    # 5. Map back to original
    df['assigned_groups'] = [[] for _ in range(len(df))]
    df['minimal_species_per_group'] = [[] for _ in range(len(df))]
    # Store unique count in a temporary dict for the collapse function
    group_unique_map = {}
    group_dfindex_map = {}

    for g_idx, group in enumerate(final_valid_groups):
        g_label = f"BACT_G_{g_idx}"
        group_unique_map[g_label] = group['n_unique_loci']
        group_dfindex_map[g_label] = group['indices']
        for row_idx in group['indices']:
            df.at[row_idx, 'assigned_groups'].append(g_label)
            df.at[row_idx, 'minimal_species_per_group'].append(group['common'])

    return df, group_unique_map,group_dfindex_map


def replace_with_collapsed_groups(df, group_unique_map, sum_cols,carryovercols,group_dfindex_map,rmlst_stats):

    if 'assigned_groups' not in df.columns:
        df = df.copy()
        df['assigned_groups'] = [[] for _ in range(len(df))]

    def _has_assigned(x):
        if isinstance(x, (list, tuple, set)):
            return len(x) > 0
        return False

    mask_used = df['assigned_groups'].apply(_has_assigned)
    exploded_df = df[mask_used].explode('assigned_groups')
    exploded_df['locus'] = exploded_df["graph_name"].str.split("_").str[0].str.replace("bac","BACT")
    # Aggregation
    group_to_species = {}
    for idx, row in df[mask_used].iterrows():
        for i, g_id in enumerate(row['assigned_groups']):
            if g_id not in group_to_species:
                group_to_species[g_id] = row['minimal_species_per_group'][i]
    # for each assigned group check for any missing loci in rmlst_stats, if they are missing add a dummy row with 0 in all sum_cols except npos_max_probetype which should be rmlst_stats[locus]["median"]

    rmlstorder = list(sorted(rmlst_stats.keys()))

    for g_id, species in group_to_species.items():
        explodedgroup = exploded_df.loc[exploded_df['assigned_groups'] == g_id]
        includedloci = explodedgroup['locus'].apply(lambda x: x in rmlstorder)
        for locus in rmlstorder:
            if locus not in explodedgroup['locus'].values:
                new_row = pd.Series({col: 0 for col in sum_cols})
                new_row['assigned_groups'] = g_id
                new_row['locus'] = locus
                if "npos_max_probetype" in sum_cols:
                    new_row["npos_max_probetype"] = rmlst_stats[locus]["median"]
                exploded_df = pd.concat([exploded_df, new_row.to_frame().T], ignore_index=True)
        ...
    for g_id, species in group_to_species.items():
        df_indices = group_dfindex_map[g_id]
        for idx in df_indices:
            row = df.loc[idx]
            if row["graph_name"] in rmlst_stats:
                locus = row["graph_name"]
                for col in sum_cols:
                    if col not in row or pd.isna(row[col]):
                        if col == "npos_max_probetype":
                            df.at[idx, col] = rmlst_stats[locus]["median"]
                        else:
                            df.at[idx, col] = 0

    agg_logic = {col: 'sum' for col in sum_cols if col in df.columns}
    for x in carryovercols: agg_logic[x] = 'first'
    agg_logic['graph_name'] = 'count'
    # Include pathlen in sum_cols if it isn't already there for the prop calculation
    if 'pathlen' not in agg_logic: agg_logic['pathlen'] = 'sum'

    collapsed_df = exploded_df.groupby('assigned_groups').agg(agg_logic).reset_index()

    # Metadata and Metrics
    collapsed_df = collapsed_df.rename(columns={'assigned_groups': 'group_id', 'graph_name': 'n_loci'})
    collapsed_df['species_set'] = collapsed_df['group_id'].map(group_to_species)
    collapsed_df['n_unique_loci'] = collapsed_df['group_id'].map(group_unique_map)
    collapsed_df["maindfindex"] = collapsed_df['group_id'].map(group_dfindex_map)
    for i in [1, 2, 5, 10, 100, 1000]:
        col_name = f"npos_cov_mindepth{i}"
        if col_name in collapsed_df.columns:
            collapsed_df[f"prop_npos_cov{i}"] = collapsed_df[col_name] / collapsed_df["npos_max_probetype"]

    # Combine
    df_remaining = df[~mask_used].copy()
    final_df = pd.concat([df_remaining, collapsed_df], ignore_index=True)


    return collapsed_df

def merge_nodes(topgraphdata,depthfiles,consensusfiles,plotdir,consdir,sampleid,rmlst = False):
    ...
    for i in topgraphdata.iterrows():
        graphname = i[1]["graph_name"]
        speciesset = i[1]["species_set"]
        speciessetls = speciesset.split(", ")
        if len(speciessetls) > 3:
            genera = list(set([x.split("~")[0] + f"-{speciesset.count(x.split("~")[0])}" for x in speciessetls]))
            genera = ", ".join(genera)
            ...
        else:
            genera = speciesset
        outheader = f"{i[1]["graph_name"].replace("_graph", "")}_{genera}"
        blocktuple = list(map(str, i[1]["block_id_order"]))
        blockcounttuple = list(i[1]["blockcounts"])
        blocklentuple = list(i[1]["block_len_tuple"])
        if not graphname.startswith("bac000"):
            blockstart = 0
            pathoverall = None
            consensusseq = ""
            for blockindex, block in enumerate(blocktuple):
                parsedblock = block.split("block")[-1].split("-")[0]
                if parsedblock in depthfiles:
                    blockdepths = pd.read_csv(depthfiles[parsedblock], header=None)
                    blockdepths.index = ["all_reads", "dedup_reads"]
                    blockdepths = blockdepths.transpose()
                    blockdepths = blockdepths.iloc[:-1]
                    blockdepths["all_reads_per_repitition"] = blockdepths["all_reads"] / blockcounttuple[blockindex]
                    blockdepths["dedup_reads_per_repitition"] = blockdepths["dedup_reads"] / blockcounttuple[blockindex]
                    blockdepths["position"] = blockdepths.index.astype(int)
                    blockdepths["pathposition"] = blockdepths["position"] + blockstart
                    blockstart = blockdepths["pathposition"].max() + 1
                    blockdepths["blockmissing"] = False
                    if pathoverall is None:
                        pathoverall = blockdepths
                    else:
                        pathoverall = pd.concat([pathoverall, blockdepths], ignore_index=False)

                else:
                    blockdepths = pd.DataFrame(columns=["position", "all_reads", "all_reads_per_repitition", "dedup_reads",
                                                        "dedup_reads_per_repitition", "pathposition", "blockmissing"])
                    blockdepths["position"] = range(blocklentuple[blockindex])
                    blockdepths["pathposition"] = blockdepths["position"] + blockstart
                    blockstart = blockdepths["pathposition"].max() + 1
                    blockdepths["all_reads"] = 0
                    blockdepths["dedup_reads"] = 0
                    blockdepths["dedup_reads_per_repitition"] = 0
                    blockdepths["all_reads_per_repitition"] = 0
                    blockdepths["blockmissing"] = True
                    if pathoverall is None:
                        pathoverall = blockdepths
                    else:
                        pathoverall = pd.concat([pathoverall, blockdepths], ignore_index=False)
                if parsedblock in consensusfiles:
                    inf = SeqIO.to_dict(SeqIO.parse(consensusfiles[parsedblock], "fasta"))
                    key = list(inf.keys())[0]
                    consensusseq += str(inf[key].seq)
                else:
                    consensusseq += "N" * blocklentuple[blockindex]

            pathoverall["pathposition"] = pd.to_numeric(pathoverall["pathposition"], errors="coerce")
            pathoverall["blockmissing"] = pathoverall["blockmissing"].astype(bool)

            subset = pathoverall[pathoverall["position"] == 0]
            x_values = subset["pathposition"].unique()
            n_all = i[1].get("n_reads_all", None)
            n_all_perrep = i[1].get("n_reads_all_per_repetition", None)
            n_dedup_perrep = i[1].get("n_reads_dedup_per_repetition", None)
            n_dedup = i[1].get("n_reads_dedup", None)
            nc2 = i[1].get("npos_cov_mindepth2", None)
            pn = i[1].get("prop_npos_cov2", None)
            pn_str = f"{pn:.3f}" if (pn is not None and pd.notna(pn)) else str(pn)
            stats_text = (
                f"n_reads_all: {n_all}\n"
                f"n_all_perrep: {n_all_perrep}\n"
                f"n_reads_dedup: {n_dedup}\n"
                f"n_dedup_perrep: {n_dedup_perrep}\n"
                f"npos_cov_mindepth2: {nc2}\n"
                f"prop_npos_cov2: {pn_str}"
            )

            if float(pn) > 0.35:

                fig, ax = plt.subplots(figsize=(30, 6))
                sns.lineplot(data=pathoverall, x="pathposition", y="dedup_reads", color="red", ax=ax)
                sns.lineplot(data=pathoverall, x="pathposition", y="dedup_reads_per_repitition", color="red", ax=ax,
                             alpha=0.3)
                sns.lineplot(data=pathoverall, x="pathposition", y="all_reads", color="blue", ax=ax)
                sns.lineplot(data=pathoverall, x="pathposition", y="all_reads_per_repitition", color="blue", ax=ax,
                             alpha=0.3)

                plt.title(f"{graphname}:{genera}")
                for x in x_values:
                    plt.axvline(x=x, linestyle="--")

                # Make room on the right and place the text
                fig.subplots_adjust(right=0.75)
                ax.text(1.02, 0.95, stats_text, transform=ax.transAxes, ha="left", va="top",
                        fontsize=10, family="monospace", bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))
                plt.yscale("log")

                ymin, ymax = ax.get_ylim()

                ax.fill_between(
                    pathoverall["pathposition"],
                    ymin,
                    ymax,
                    where=pathoverall["blockmissing"],
                    alpha=0.7,
                    color="grey"
                )

                ax.set_ylim(ymin, ymax)
                sgraphname = re.sub(r'[^A-Za-z0-9._-]', '_', graphname).replace("_graph", "")
                sspeciesset = re.sub(r'[^A-Za-z0-9._-]', '_', genera)
                filename = f"{plotdir}/{sgraphname}-{sspeciesset}_log.png"

                plt.savefig(filename, dpi=1000)
                plt.close()

                fig, ax = plt.subplots(figsize=(30, 6))
                sns.lineplot(data=pathoverall, x="pathposition", y="dedup_reads", color="red", ax=ax)
                sns.lineplot(data=pathoverall, x="pathposition", y="dedup_reads_per_repitition", color="red", ax=ax,
                             alpha=0.3)
                sns.lineplot(data=pathoverall, x="pathposition", y="all_reads", color="blue", ax=ax)
                sns.lineplot(data=pathoverall, x="pathposition", y="all_reads_per_repitition", color="blue", ax=ax,
                             alpha=0.3)

                plt.title(f"{graphname}:{genera}")
                for x in x_values:
                    plt.axvline(x=x, linestyle="--")

                # Make room on the right and place the text
                fig.subplots_adjust(right=0.75)
                ax.text(1.02, 0.95, stats_text, transform=ax.transAxes, ha="left", va="top",
                        fontsize=10, family="monospace", bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))
                # plt.yscale("log")

                ymin, ymax = ax.get_ylim()

                ax.fill_between(
                    pathoverall["pathposition"],
                    ymin,
                    ymax,
                    where=pathoverall["blockmissing"],
                    alpha=0.7,
                    color="grey"
                )

                ax.set_ylim(ymin, ymax)
                sgraphname = re.sub(r'[^A-Za-z0-9._-]', '_', graphname).replace("_graph", "")
                sspeciesset = re.sub(r'[^A-Za-z0-9._-]', '_', genera)
                filename = f"{plotdir}/{sgraphname}-{sspeciesset}.png"

                plt.savefig(filename, dpi=1000)
                plt.close()

                nuccount = sum(consensusseq.count(b) for b in "ATGCatgc")
                if nuccount > 0:
                    graphspecs = re.sub(r'[^A-Za-z0-9._-]', '_', graphname).replace("_graph", "") + "_" + re.sub(
                        r'[^A-Za-z0-9._-]', '_', genera)
                    pathconsensus = SeqRecord.SeqRecord(Seq.Seq(consensusseq), id=f"{graphspecs}_{sampleid}_consensus",
                                                        description="")
                    SeqIO.write(pathconsensus, f"{consdir}/{graphspecs}_{sampleid}_consensus.fasta", "fasta")


def plotrmlst(pathoverall,graphname,genera,x_values,locuspos,stats_text,plotdir,islog=False):
    fig, ax = plt.subplots(figsize=(30, 6))
    sns.lineplot(data=pathoverall, x="locusposition", y="dedup_reads", color="red", ax=ax)
    sns.lineplot(data=pathoverall, x="locusposition", y="dedup_reads_per_repitition", color="red", ax=ax,
                 alpha=0.3)
    sns.lineplot(data=pathoverall, x="locusposition", y="all_reads", color="blue", ax=ax)
    sns.lineplot(data=pathoverall, x="locusposition", y="all_reads_per_repitition", color="blue", ax=ax,
                 alpha=0.3)

    plt.title(f"{graphname}:{genera}")
    for x in x_values:
        plt.axvline(x=x, linestyle="--", alpha=0.2)
    for x in locuspos:
        plt.axvline(x=x, linestyle="--")
    # Make room on the right and place the text
    fig.subplots_adjust(right=0.75)
    ax.text(1.02, 0.95, stats_text, transform=ax.transAxes, ha="left", va="top",
            fontsize=10, family="monospace", bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))
    if islog:
        plt.yscale("log")


    ymin, ymax = ax.get_ylim()

    ax.fill_between(
        pathoverall["locusposition"],
        ymin,
        ymax,
        where=pathoverall["blockmissing"],
        alpha=0.7,
        color="grey"
    )

    ax.set_ylim(ymin, ymax)
    sgraphname = re.sub(r'[^A-Za-z0-9._-]', '_', graphname).replace("_graph", "")
    sspeciesset = re.sub(r'[^A-Za-z0-9._-]', '_', genera)
    if islog:
        filename = f"{plotdir}/{sgraphname}-{sspeciesset}_log.png"
    else:
        filename = f"{plotdir}/{sgraphname}-{sspeciesset}.png"

    plt.savefig(filename, dpi=1000)
    plt.close()

def merge_rmlst_nodes(rmlstgroups,rmlstgraphres,depthfiles,consensusfiles,outdir,sampleid):
    ...
    #Similar to merge_nodes but for each rmlst locus in a group generate a folder where each locus in that group has consensus and plots
    # then generate a combined plot across all loci also generate a combined fasta file
    rmlstdir = f"{outdir}/rmlst/"
    plotdir = f"{rmlstdir}/plots/"
    consdir = f"{rmlstdir}/consensus/"
    if not os.path.isdir(rmlstdir):
        os.makedirs(rmlstdir)
        os.mkdir(consdir)
        os.mkdir(plotdir)
    else:
        shutil.rmtree(rmlstdir)
        os.makedirs(rmlstdir)
        os.mkdir(consdir)
        os.mkdir(plotdir)

    for i in rmlstgroups.iterrows():
        indices = i[1]["maindfindex"]

        speciesset = i[1]["species_set"]
        speciessetls = list(speciesset)
        if len(speciessetls) > 3:
            genera = list(set([x.split("~")[0] + f"-{speciesset.count(x.split("~")[0])}" for x in speciessetls]))
            genera = "-".join(genera)
            ...
        else:
            genera = speciesset
            genera = "-".join(genera)
        outheader = f"{i[1]["graph_name"].replace("_graph", "")}_{genera}"
        # groupdir = f"{rmlstdir}/{outheader}"
        # os.mkdir(groupdir)
        # plotdir = f"{groupdir}/plots/"
        # consdir = f"{groupdir}/consensus/"

        rmlstorder = list(sorted(rmlst_stats.keys()))

        locus_to_indices = {rmlstgraphres.iloc[index]["graph_name"].replace("bac","BACT").replace("_graph",""):index for index in indices}
        indicesorder = [locus_to_indices[x] if x in locus_to_indices.keys() else None for x in rmlstorder ]
        pathoverall = None
        consensusseq = ""
        locusstart = 0
        consensusseqs = []
        for locus in rmlstorder:
            blockstart = 0
            mediansize = rmlst_stats[locus]["median"]
            if locus in locus_to_indices:
                # print(locus,"graph")
                index = locus_to_indices[locus]
                grouplocus = rmlstgraphres.iloc[index]
                graphname = grouplocus["graph_name"]
                graphname = graphname.split("_")[0]
                if grouplocus["isgraph"]:
                # graphname = grouplocus["graph_name"]
                    blocktuple = list(map(str, grouplocus["block_id_order"]))
                    blockcounttuple = list(grouplocus["blockcounts"])
                    blocklentuple = list(grouplocus["block_len_tuple"])
                    # if not graphname.startswith("bac000"):

                    for blockindex, block in enumerate(blocktuple):
                        parsedblock = block.split("block")[-1].split("-")[0]
                        if parsedblock in depthfiles:
                            blockdepths = pd.read_csv(depthfiles[parsedblock], header=None)
                            blockdepths.index = ["all_reads", "dedup_reads"]
                            blockdepths = blockdepths.transpose()
                            blockdepths = blockdepths.iloc[:-1]
                            blockdepths["all_reads_per_repitition"] = blockdepths["all_reads"] / blockcounttuple[blockindex]
                            blockdepths["dedup_reads_per_repitition"] = blockdepths["dedup_reads"] / blockcounttuple[blockindex]
                            blockdepths["position"] = blockdepths.index.astype(int)
                            blockdepths["locusposition"] = blockdepths["position"] + locusstart
                            blockdepths["pathposition"] = blockdepths["position"] + blockstart
                            blockstart = blockdepths["pathposition"].max() + 1
                            locusstart = blockdepths["locusposition"].max() + 1
                            blockdepths["blockmissing"] = False
                            if pathoverall is None:
                                pathoverall = blockdepths
                            else:
                                pathoverall = pd.concat([pathoverall, blockdepths], ignore_index=False)

                        else:
                            blockdepths = pd.DataFrame(
                                columns=["position", "all_reads", "all_reads_per_repitition", "dedup_reads",
                                         "dedup_reads_per_repitition", "pathposition", "blockmissing","locusposition"])
                            blockdepths["position"] = range(blocklentuple[blockindex])
                            blockdepths["pathposition"] = blockdepths["position"] + blockstart
                            blockdepths["locusposition"] = blockdepths["position"] + locusstart
                            blockstart = blockdepths["pathposition"].max() + 1
                            locusstart = blockdepths["locusposition"].max() + 1
                            blockdepths["all_reads"] = 0
                            blockdepths["dedup_reads"] = 0
                            blockdepths["dedup_reads_per_repitition"] = 0
                            blockdepths["all_reads_per_repitition"] = 0
                            blockdepths["blockmissing"] = True
                            if pathoverall is None:
                                pathoverall = blockdepths
                            else:
                                pathoverall = pd.concat([pathoverall, blockdepths], ignore_index=False)
                        if parsedblock in consensusfiles:
                            inf = SeqIO.to_dict(SeqIO.parse(consensusfiles[parsedblock], "fasta"))
                            key = list(inf.keys())[0]
                            consensusseq += str(inf[key].seq)
                        else:
                            consensusseq += "N" * blocklentuple[blockindex]
                else:
                    # print(locus, "nongraph")
                    #TODO resolve below, need to use single node as full path
                    blockid = grouplocus["block_idonly"]
                    blockdepths = pd.read_csv(depthfiles[blockid], header=None)
                    blockdepths.index = ["all_reads", "dedup_reads"]
                    blockdepths = blockdepths.transpose()
                    blockdepths = blockdepths.iloc[:-1]
                    blockdepths["all_reads_per_repitition"] = blockdepths["all_reads"] / blockcounttuple[blockindex]
                    blockdepths["dedup_reads_per_repitition"] = blockdepths["dedup_reads"] / blockcounttuple[blockindex]
                    blockdepths["position"] = blockdepths.index.astype(int)
                    blockdepths["locusposition"] = blockdepths["position"] + locusstart
                    blockdepths["pathposition"] = blockdepths["position"] + blockstart
                    blockstart = blockdepths["pathposition"].max() + 1
                    locusstart = blockdepths["locusposition"].max() + 1
                    blockdepths["blockmissing"] = False
                    if pathoverall is None:
                        pathoverall = blockdepths
                    else:
                        pathoverall = pd.concat([pathoverall, blockdepths], ignore_index=False)

                    # blockdepths = pd.DataFrame(
                    #     columns=["position", "all_reads", "all_reads_per_repitition", "dedup_reads",
                    #              "dedup_reads_per_repitition", "pathposition", "blockmissing"])
                    # blockdepths["position"] = range(mediansize)
                    # blockdepths["pathposition"] = blockdepths["position"] + blockstart
                    # blockdepths["locusposition"] = blockdepths["position"] + locusstart
                    # blockstart = blockdepths["pathposition"].max() + 1
                    # locusstart = blockdepths["locusposition"].max() + 1
                    # blockdepths["all_reads"] = 0
                    # blockdepths["dedup_reads"] = 0
                    # blockdepths["dedup_reads_per_repitition"] = 0
                    # blockdepths["all_reads_per_repitition"] = 0
                    # blockdepths["blockmissing"] = True
                    # if pathoverall is None:
                    #     pathoverall = blockdepths
                    # else:
                    #     pathoverall = pd.concat([pathoverall, blockdepths], ignore_index=False)

            else:
                # print(locus,"missing")
                blockdepths = pd.DataFrame(
                    columns=["position", "all_reads", "all_reads_per_repitition", "dedup_reads",
                             "dedup_reads_per_repitition", "pathposition", "blockmissing"])
                blockdepths["position"] = range(mediansize)
                blockdepths["pathposition"] = blockdepths["position"] + blockstart
                blockdepths["locusposition"] = blockdepths["position"] + locusstart
                blockstart = blockdepths["pathposition"].max() + 1
                locusstart = blockdepths["locusposition"].max() + 1
                blockdepths["all_reads"] = 0
                blockdepths["dedup_reads"] = 0
                blockdepths["dedup_reads_per_repitition"] = 0
                blockdepths["all_reads_per_repitition"] = 0
                blockdepths["blockmissing"] = True
                if pathoverall is None:
                    pathoverall = blockdepths
                else:
                    pathoverall = pd.concat([pathoverall, blockdepths], ignore_index=False)

            nuccount = sum(consensusseq.count(b) for b in "ATGCatgc")
            if nuccount > 0:
                graphspecs = re.sub(r'[^A-Za-z0-9._-]', '_', graphname).replace("_graph", "") + "_" + re.sub(
                    r'[^A-Za-z0-9._-]', '_', genera)
                pathconsensus = SeqRecord.SeqRecord(Seq.Seq(consensusseq), id=f"{graphspecs}_{sampleid}_consensus",
                                                    description="")
                consensusseqs.append(pathconsensus)


        pathoverall["pathposition"] = pd.to_numeric(pathoverall["pathposition"], errors="coerce")
        pathoverall["pathposition"] = pd.to_numeric(pathoverall["pathposition"], errors="coerce")
        pathoverall["blockmissing"] = pathoverall["blockmissing"].astype(bool)

        subset = pathoverall[pathoverall["position"] == 0]
        x_values = subset["pathposition"].unique()
        locussubset = pathoverall[pathoverall["pathposition"] == 0]
        locuspos = locussubset["locusposition"].unique()

        groupid = i[1]["group_id"]
        grouplocus = rmlstgroups.loc[rmlstgroups["group_id"] == groupid]
        n_all = grouplocus["n_reads_all"].tolist()[0]
        # n_all_perrep = grouplocus.get("n_reads_all_per_repetition", None)
        # n_dedup_perrep = grouplocus.get("n_reads_dedup_per_repetition", None)
        n_dedup = grouplocus["n_reads_dedup"].tolist()[0]
        nc2 = grouplocus["npos_cov_mindepth2"].tolist()[0]
        pn = grouplocus["prop_npos_cov2"].tolist()[0]

        pn_str = f"{pn:.3f}" if (pn is not None and pd.notna(pn)) else str(pn)
        stats_text = (
            f"n_reads_all: {n_all:0.0f}\n"
            f"n_reads_dedup: {n_dedup:0.0f}\n"
            f"npos_cov_mindepth2: {nc2:0.0f}\n"
            f"prop_npos_cov2: {pn_str}"
        )
        ...
        if float(pn) > 0.35:

            plotrmlst(pathoverall,"rMLST", genera, x_values, locuspos, stats_text, plotdir)
            plotrmlst(pathoverall, "rMLST", genera, x_values, locuspos, stats_text, plotdir,islog=True)



            #TODO rather than doing plots for each locus, generate one plots across all rMLST loci
            #TODO also generate one consensus fasta for the sample/rmlst group
            ...
        SeqIO.write(consensusseqs, f"{consdir}/{sampleid}_{outheader}_consensus.fasta", "fasta")
    ...




def main():
    """
    Main processing pipeline:
    1. Read input graph dataframe (pickle) and sample Castanet CSV.
    2. Derive a 'block_id' for each probe row by removing the 'rmlst' prefix.
    3. Map per-block coverage/read metrics onto each path (sum over the block set).
    4. Compute proportions and filter paths by read support and coverage proportion.
    5. Remove paths that are strict subsets of a higher-ranked path (by support/proportion).
    6. Select top path per (component, graph_name).
    """
    args = get_args()

    # Load inputs provided by the user via CLI
    # graphdata = pd.read_pickle(, compression="gzip")
    graphdata = pd.read_parquet(args.inputgraphdata)

    graphdata["block_id_set"] = graphdata.apply(update_to_longblocksetnames,axis=1)
    graphdata["block_id_order"] = graphdata.apply(update_to_longblockordernames,axis=1)
    graphdata["graph_name"] = graphdata["graph_name"].str.replace("BACT000","bac000").replace("16S","16s").replace("23S","23s")
    graphdata["blockcounts"] = graphdata["block_id_order"].apply(lambda lst: [Counter(lst)[x] for x in lst])
    species_dict = {
        s.replace("~", "").lower(): s
        for sublist in graphdata["species_set"]
        for s in sublist
    }

    samplecastanet = pd.read_csv(args.inputcastanet)

    if "n_reads_all" not in samplecastanet.columns:
        samplecastanet["n_reads_all"] = samplecastanet["reads_for_mapping"]
        samplecastanet["clean_n_reads_all"] = samplecastanet["clean_reads_for_mapping"]

    sampleid = samplecastanet['sampleid'].tolist()[0]
    # Quick sanity check: ensure expected columns exist in the Castanet CSV
    expected_cols = {'probetype', 'npos_cov_mindepth2', 'n_reads_dedup'}
    missing = expected_cols.difference(set(samplecastanet.columns))
    if missing:
        raise ValueError(f"Input castanet file is missing required columns: {missing}")
    all_block_ids = []



    all_block_ids = (
        graphdata["block_id_set"]
        .dropna()
        .explode()
        .unique()
        .tolist()
    )

    samplecastanet["block_id"] = samplecastanet["probetype"]
    samplecastanet["block_idonly"] = samplecastanet["probetype"].str.split("block").str[-1].str.split("_").str[0].str.split("-").str[0]
    samplecastanet["graph_id"] = samplecastanet["probetype"].str.split("block").str[0].replace("graph", "")

    #TODO when matching blocks from castanet and graphdata need to match graph and block - some blocks have same block_id but are in different graphs

    non_graph_hits = samplecastanet[~samplecastanet["block_id"].isin(all_block_ids)].copy()

    #TODO only use node a present if it is over 0.35 covered
    sumcols = ["npos_cov_mindepth2", "n_reads_dedup","n_reads_all",'npos_max_probetype',
       'npos_cov_probetype','npos_cov_mindepth1', 'npos_cov_mindepth2', 'npos_cov_mindepth5',
       'npos_cov_mindepth10', 'npos_cov_mindepth100', 'npos_cov_mindepth1000','npos_dedup_cov_mindepth1', 'npos_dedup_cov_mindepth2',
       'npos_dedup_cov_mindepth5', 'npos_dedup_cov_mindepth10',
       'npos_dedup_cov_mindepth100', 'npos_dedup_cov_mindepth1000']
    sumcolswiththresh = ["npos_cov_mindepth2", "n_reads_dedup","n_reads_all",'npos_max_probetype',
       'npos_cov_probetype','npos_cov_mindepth1', 'npos_cov_mindepth2', 'npos_cov_mindepth5',
       'npos_cov_mindepth10', 'npos_cov_mindepth100', 'npos_cov_mindepth1000','npos_dedup_cov_mindepth1', 'npos_dedup_cov_mindepth2',
       'npos_dedup_cov_mindepth5', 'npos_dedup_cov_mindepth10',
       'npos_dedup_cov_mindepth100', 'npos_dedup_cov_mindepth1000']
    divbycoount = ["n_reads_dedup","n_reads_all"]
    carryovercols = ['reads_on_target', 'reads_on_target_dedup',]

    for col in carryovercols:
        if col not in samplecastanet.columns:
            raise ValueError(f"Input castanet file is missing required column: {col}")
        colmap = samplecastanet.set_index("block_id")[col].to_dict()
        graphdata[col] = graphdata["block_id_order"].str[0]
    for col in sumcols:
        if col not in samplecastanet.columns:
            raise ValueError(f"Input castanet file is missing required column: {col}")
        colmap = samplecastanet.set_index("block_id")[col].to_dict()
        graphdata[col] = graphdata["block_id_order"].apply(
        lambda s: sum(colmap.get(bid, 0) for bid in s)
        )
        if col in divbycoount:
            modcol = col+"_per_repetition"
            #TODO fix below

            graphdata[modcol] = graphdata.apply(lambda row: row[col] / sum(row["blockcounts"]) if sum(row["blockcounts"]) > 0 else 0, axis=1)
    for col in sumcolswiththresh:
        if col not in samplecastanet.columns:
            raise ValueError(f"Input castanet file is missing required column: {col}")
        colmap = samplecastanet.set_index("block_id")[col].to_dict()
        threshmap = samplecastanet.set_index("block_id")[f"prop_npos_cov2"].to_dict()
        graphdata[col+"_thresh"] = graphdata["block_id_order"].apply(
        lambda s: sum(colmap.get(bid, 0) for bid in s if threshmap.get(bid, 0) > 0.35)
        )
        if col in divbycoount:
            modcol = col+"_per_repetition"
            graphdata[modcol] = graphdata.apply(
                lambda row: sum(
                    colmap.get(bid, 0) / cnt for bid, cnt in zip(row["block_id_order"], row["blockcounts"])) if sum(
                    row["blockcounts"]) > 0 else 0,
                axis=1,
            )
            # graphdata[modcol] = graphdata.apply(lambda row: row[col] / sum(row["blockcounts"]) if sum(row["blockcounts"]) > 0 else 0, axis=1)

    meancols = ['amprate_mean','depth_mean','udepth_mean']
    for col in meancols:
        if col not in samplecastanet.columns:
            raise ValueError(f"Input castanet file is missing required column: {col}")
        colmap = samplecastanet.set_index("block_id")[col].to_dict()
        graphdata[col] = graphdata["block_id_order"].apply(
        lambda s: statistics.mean([colmap.get(bid, 0) for bid in s])
        )

    usethresh = False
    # Compute the proportion of positions covered at depth>=2 relative to path length
    if usethresh:
        for i in [1,2,5,10,100,1000]:
            graphdata[f"prop_npos_cov{i}"] = graphdata[f"npos_cov_mindepth{i}_thresh"] / graphdata["npos_max_probetype"]
    else:
        for i in [1, 2, 5, 10, 100, 1000]:
            graphdata[f"prop_npos_cov{i}"] = graphdata[f"npos_cov_mindepth{i}"] / graphdata["npos_max_probetype"]

    # Filter out paths with no dedup reads and low coverage proportion
    graphdata = graphdata[graphdata["n_reads_dedup"] > 0]
    # graphdata = graphdata[graphdata["prop_npos_cov2"] > 0.1]

    #if no rows remain after filtering, exit
    if graphdata.shape[0]==0:
        print("No paths remain after filtering by read support and coverage proportion. Exiting.")
        outdf = pd.DataFrame(columns=["graph_name","component","pathlen","npos_cov_mindepth2","cumulative_dedup_reads","prop_npos_cov2","block_id_order","species_set"])
        outdf.to_csv(f"{args.output}_top_paths.tsv",sep="\t", index=False)
        return
    # For each graph_name, remove any path that is a strict subset of a higher-ranked path
    # Ranking is by prop_npos_cov2 (both descending)
    finalgraphdata = pd.DataFrame()
    for gname, gdf in graphdata.groupby("graph_name"):
        # Sort so that higher-support / higher-proportion paths appear earlier
        gdf = gdf.sort_values(by=["prop_npos_cov2"], ascending=False).reset_index(drop=True)
        todrop = set()
        # Compare each path to later ones: if current path's block set is a superset of a later
        # path's block set, then the later path is considered a subset and dropped (it has less support)
        for i in range(len(gdf)):
            if i in todrop:
                # already marked for removal
                continue
            bi = gdf.loc[i, "block_id_set"]
            bi = set(bi)
            for j in range(i + 1, len(gdf)):
                if j in todrop:
                    continue
                bj = gdf.loc[j, "block_id_set"]
                bj = set(bj)
                # If the earlier (better-ranked) path covers all blocks of a later path,
                # then drop the later path as redundant
                if bi.issuperset(bj):
                    todrop.add(j)
        finalgdf = gdf.drop(index=todrop).reset_index(drop=True)
        finalgraphdata = pd.concat([finalgraphdata, finalgdf], axis=0)

    # From the cleaned set, choose the top path per (component, graph_name) by the same ranking
    finalgraphdata = finalgraphdata.sort_values(by=["npos_cov_mindepth2"], ascending=False).reset_index(
        drop=True)
    nonrmlstfinal = finalgraphdata[~finalgraphdata['graph_name'].str.startswith('bac000', na=False)]
    # Restore grouping by component and graph_name to select the top (highest-ranked) path
    #TODO forS87 in 16S wrong path being selected.
    check = finalgraphdata.groupby(["component", "graph_name"], sort=False).count().reset_index()
    topgraphdata = finalgraphdata.groupby(["component", "graph_name"], sort=False).first().reset_index()

    #process rMLST data

    mask = topgraphdata['graph_name'].str.startswith('bac000', na=False)
    rmlstdata = topgraphdata[mask].copy()

    nongraph_rmlstdata = non_graph_hits[non_graph_hits["graph_id"].str.startswith('graphbac000', na=False)].copy()
    # nongraph_rmlstdata["graph_name"] = nongraph_rmlstdata["graph_id"].str.replace("graph", "")
    nongraph_rmlstdata[["graph_name", "species_set"]] = nongraph_rmlstdata["probetype"].str.extract(r"block(bac\d+)(.+)$")
    nongraph_rmlstdata["species_set"] = nongraph_rmlstdata["species_set"].map(lambda x: [species_dict[x]] if x in species_dict else [x])
    # nongraph_rmlstdata["species_set"] = nongraph_rmlstdata["species_set"].map(lambda x: [x])
    nongraph_rmlstdata["isgraph"] = False
    rmlstdata["isgraph"] = True
    merged_rmlst = pd.concat([rmlstdata, nongraph_rmlstdata], ignore_index=True).reset_index()
    rmlstdata,group_metadata,group_dfindex_map = process_strict_unique_bact_groups(merged_rmlst)

    if group_metadata != {}:
        rmlstsumcols = ["pathlen","npos_cov_mindepth2", "n_reads_dedup","n_reads_all",'npos_max_probetype',
           'npos_cov_probetype','npos_cov_mindepth1', 'npos_cov_mindepth2', 'npos_cov_mindepth5',
           'npos_cov_mindepth10', 'npos_cov_mindepth100', 'npos_cov_mindepth1000','npos_dedup_cov_mindepth1', 'npos_dedup_cov_mindepth2',
           'npos_dedup_cov_mindepth5', 'npos_dedup_cov_mindepth10',
           'npos_dedup_cov_mindepth100', 'npos_dedup_cov_mindepth1000',]
        carryovercols = ['block_id_order',"block_len_tuple",""'reads_on_target', 'reads_on_target_dedup',]
        rmlstdata = replace_with_collapsed_groups(merged_rmlst,group_metadata,rmlstsumcols,carryovercols,group_dfindex_map,rmlst_stats)
        rmlstdata["graph_name"] = "rmlst_"+rmlstdata["group_id"].str.replace("BACT_","")
        #group_id species_set n_loci n_unique_loci pathlen	npos_cov_mindepth2	n_reads_dedup	n_reads_all	npos_max_probetype	npos_cov_probetype	npos_cov_mindepth1	npos_cov_mindepth5	npos_cov_mindepth10	npos_cov_mindepth100	npos_cov_mindepth1000	npos_dedup_cov_mindepth1	npos_dedup_cov_mindepth2	npos_dedup_cov_mindepth5	npos_dedup_cov_mindepth10	npos_dedup_cov_mindepth100	npos_dedup_cov_mindepth1000 prop_npos_cov1	prop_npos_cov2	prop_npos_cov5	prop_npos_cov10	prop_npos_cov100	prop_npos_cov1000
        colforfinal = ['group_id','species_set','n_loci','n_unique_loci','pathlen','npos_cov_mindepth2','n_reads_dedup','n_reads_all','npos_max_probetype','npos_cov_probetype','npos_cov_mindepth1','npos_cov_mindepth5','npos_cov_mindepth10','npos_cov_mindepth100','npos_cov_mindepth1000','npos_dedup_cov_mindepth1','npos_dedup_cov_mindepth2','npos_dedup_cov_mindepth5','npos_dedup_cov_mindepth10','npos_dedup_cov_mindepth100','npos_dedup_cov_mindepth1000', 'prop_npos_cov1', 'prop_npos_cov2', 'prop_npos_cov5', 'prop_npos_cov10', 'prop_npos_cov100', 'prop_npos_cov1000']
        rmlstdatafinal = rmlstdata[colforfinal]
        rmlstdatafinal["species_set"] = rmlstdatafinal["species_set"].apply(lambda x: ", ".join(x))
        rmlstdatafinal.to_csv(f"{args.output}_rMLST_combined.tsv",sep="\t", index=False)
        rmlstdata.to_csv(f"{args.output}_rMLST_debug.tsv", sep="\t", index=False)




    # topgraphdata = topgraphdata[~mask]

    topgraphdata["species_set"] = topgraphdata.apply(
        update_species_set,
        axis=1
    )

    non_graph_hits = non_graph_hits.copy()
    non_graph_hits["graph_name"] = non_graph_hits["probetype"].str.split("block").str[0].str.replace("graph", "")
    non_graph_hits["species_set"] = non_graph_hits["probetype"].str.split("block").str[-1]


    topgraphdata.drop('pathnames', axis=1)
    topgraphdata = topgraphdata[~topgraphdata["graph_name"].str.fullmatch(r"\d+")]

    topgraphdata.to_csv(f"{args.output}_top_graph_paths.tsv",sep="\t", index=False)

    df_merged = topgraphdata.combine_first(non_graph_hits)
    """amprate_mean	amprate_median	amprate_std	block_id_set	clean_n_reads_all	clean_prop_of_reads_on_target	component	depth_25pc	depth_75pc	depth_mean	depth_median	depth_std	graph_name	log10_depthmean	log10_udepthmean	n_genes	n_reads_all	n_reads_dedup	n_targets	nmax_genes	nmax_targets	npos_cov_mindepth1	npos_cov_mindepth10	npos_cov_mindepth100	npos_cov_mindepth1000	npos_cov_mindepth2	npos_cov_mindepth5	npos_cov_probetype	npos_dedup_cov_mindepth1	npos_dedup_cov_mindepth10	npos_dedup_cov_mindepth100	npos_dedup_cov_mindepth1000	npos_dedup_cov_mindepth2	npos_dedup_cov_mindepth5	npos_max_probetype	path_ids	pathlen	pathnames	probetype	prop_ngenes	prop_npos_cov1	prop_npos_cov10	prop_npos_cov100	prop_npos_cov1000	prop_npos_cov2	prop_npos_cov5	prop_ntargets	prop_of_reads_on_target	pt	rawreadnum	readprop	reads_on_target	reads_on_target_dedup	sampleid	species_set	udepth_25pc	udepth_75pc	udepth_mean	udepth_median	udepth_std"""
    cols = ['graph_name', 'sampleid','component','pathlen','n_reads_dedup','n_reads_dedup_per_repetition','block_id_set','block_id_order','blockcounts','species_set',
    'probetype', 'npos_cov_mindepth1','npos_cov_mindepth2','npos_cov_mindepth5', 'npos_cov_mindepth10', 'npos_cov_mindepth100', 'npos_cov_mindepth1000',
    'npos_cov_probetype', 'npos_dedup_cov_mindepth1', 'npos_dedup_cov_mindepth2', 'npos_dedup_cov_mindepth5', 'npos_dedup_cov_mindepth10', 'npos_dedup_cov_mindepth100',
    'npos_dedup_cov_mindepth1000', 'npos_max_probetype',
    'prop_ngenes', 'prop_npos_cov1', 'prop_npos_cov2', 'prop_npos_cov5', 'prop_npos_cov10', 'prop_npos_cov100',
    'prop_npos_cov1000', 'prop_ntargets', 'prop_of_reads_on_target',
    'reads_on_target', 'reads_on_target_dedup', 'amprate_mean', 'amprate_median', 'amprate_std',
    'depth_25pc', 'depth_75pc', 'depth_mean', 'depth_median', 'depth_std',
    'log10_depthmean', 'log10_udepthmean',
    'n_genes', 'n_reads_all','n_reads_all_per_repetition', 'n_targets', 'nmax_genes', 'nmax_targets',
    'udepth_25pc', 'udepth_75pc', 'udepth_mean', 'udepth_median', 'udepth_std',
    'clean_n_reads_all', 'clean_prop_of_reads_on_target', 'pt', 'rawreadnum', 'readprop']
    df_merged = df_merged[cols]
    df_merged["sampleid"] = df_merged["sampleid"].fillna(samplecastanet["sampleid"].iloc[0])
    df_merged.to_csv(f"{args.output}_all_targets_w_top_paths.tsv", sep="\t", index=False)
    plotdir = f"{args.output}_depthplots"
    consdir = f"{args.output}_consensus"
    if not os.path.exists(plotdir):
        os.mkdir(plotdir)
    if not os.path.exists(consdir):
        os.mkdir(consdir)

    if os.path.exists(args.consensus):
        consensusfiles = {}
        for i in os.listdir(args.consensus):
            #get each folder and check that it contains a file ending in _remapped_consensus_sequence.fasta
            if os.path.isdir(os.path.join(args.consensus, i)):
                consensusfile = glob.glob(os.path.join(args.consensus, i, "*_remapped_consensus_sequence.fasta"))
                if len(consensusfile) > 0:
                    consensusfiles[i] = consensusfile[0]
        consensusfiles = {block.split("block")[-1].split("-")[0]: consensusfiles[block] for block in consensusfiles}

    if os.path.exists(args.inputdepthfolder):
        depthfiles = glob.glob(f"{args.inputdepthfolder}/*depth_by_pos.csv")
        depthfiles = {f.split("/")[-1].split("-")[0]:f for f in depthfiles}
        depthfiles = {block.split("block")[-1].split("-")[0]:depthfiles[block] for block in depthfiles}

        if group_metadata != {}:
            merge_rmlst_nodes(rmlstdata,merged_rmlst, depthfiles, consensusfiles, args.output, sampleid)


        for i in non_graph_hits.iterrows():
            if i[1].get("prop_npos_cov2", 0) > 0.1:
                # print(i[1].get("block_id", 0))
                if not i[1].get("block_id", 0).startswith("graphbac000"):
                    graphname = i[1]["probetype"]
                    # sampleid = i[1]["sampleid"]
                    outheader = f"{i[1]["graph_name"].replace("graph","")}_{i[1]['species_set']}"
                    if os.path.exists(f"{args.inputdepthfolder}/{graphname}-{sampleid}_depth_by_pos.csv"):
                        n_all = i[1].get("n_reads_all", None)

                        n_dedup = i[1].get("n_reads_dedup", None)
                        nc2 = i[1].get("npos_cov_mindepth2", None)
                        pn = i[1].get("prop_npos_cov2", None)
                        pn_str = f"{pn:.3f}" if (pn is not None and pd.notna(pn)) else str(pn)
                        stats_text = (
                            f"n_reads_all: {n_all}\n"
                            f"n_reads_dedup: {n_dedup}\n"
                            f"npos_cov_mindepth2: {nc2}\n"
                            f"prop_npos_cov2: {pn_str}"
                        )

                        if float(pn) > 0.35:
                            blockdepths = pd.read_csv(f"{args.inputdepthfolder}/{graphname}-{sampleid}_depth_by_pos.csv", header=None)
                            blockdepths.index = ["all_reads", "dedup_reads"]
                            blockdepths = blockdepths.transpose()
                            blockdepths = blockdepths.iloc[:-1]
                            blockdepths["position"] = blockdepths.index.astype(int)
                            fig, ax = plt.subplots(figsize=(30, 6))
                            sns.lineplot(data=blockdepths, x="position", y="dedup_reads", color="red", ax=ax)
                            sns.lineplot(data=blockdepths, x="position", y="all_reads", color="blue", ax=ax)
                            plt.title(f"{outheader}:{sampleid}")
                            fig.subplots_adjust(right=0.75)
                            ax.text(1.02, 0.95, stats_text, transform=ax.transAxes, ha="left", va="top",
                                    fontsize=10, family="monospace",
                                    bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))
                            plt.yscale("log")
                            filename = f"{plotdir}/{outheader}-{sampleid}_non_graph_hit.png"
                            plt.savefig(filename, dpi=1000)
                            plt.close()
                    graphcons = f"{args.consensus}/{graphname}/{graphname}_remapped_consensus_sequence.fasta"
                    if os.path.exists(graphcons):
                        s = SeqIO.parse(graphcons,"fasta")
                        for record in s:
                            consensusseq = str(record.seq)
                            nuccount = sum(consensusseq.count(b) for b in "ATGCatgc")
                            if nuccount > 0:
                                # graphspecs = re.sub(r'[^A-Za-z0-9._-]', '_', str(graphname)).replace("_graph", "") + "_" + re.sub(r'[^A-Za-z0-9._-]', '_', str(sampleid[:30]))
                                pathconsensus = SeqRecord.SeqRecord(Seq.Seq(consensusseq), id=f"{outheader}_consensus",
                                                                    description="")
                                SeqIO.write(pathconsensus, f"{consdir}/{outheader}_consensus.fasta", "fasta")
                        #,f"{consdir}/{graphname}_consensus.fasta")

        merge_nodes(topgraphdata, depthfiles, consensusfiles, plotdir, consdir, sampleid)



if __name__ == "__main__":
    main()
