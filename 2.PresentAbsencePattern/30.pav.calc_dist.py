# calculate gene distance based on gain/loss pattern
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
import umap
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA

target = ["PRKAB2", "PRDM16", "DIO2", "UCP3", "PPARGC1A", "CIDEA", "UCP1", "ATP2A1", "ACTN3"]
TM = "TMEM41B"

#%% clean data
# Read data

pav_all = pd.read_csv("data/homology/homology.pav.csv")
# gene_mm = pd.read_csv("tables/gene_DEG_top10p.csv")
# hs2mm = pd.read_csv("data/gene_name.hs2mm.csv")
# hs2mm = hs2mm.rename(columns={"gene.hs":"gene_hs", "gene.mm":"gene_mm"})
# gene_join = gene_mm.rename(columns={"gene": "gene_mm"}).merge(hs2mm, on="gene_mm", how="left")
gene_hs = pd.read_csv("tables/gene_DEG_top10p.csv")
gene_join = gene_hs.rename(columns={'gene':'gene_hs'})

pav_all = pav_all.rename(columns={"gene.hs":"gene_hs"})
gene_join = gene_join.merge(pav_all[["gene_hs"]].assign(pav_data=1), on="gene_hs",  how="left")

genes = gene_join[gene_join["pav_data"].notna()]["gene_hs"].tolist()
genes_list = list(set(genes + target))
pav = pav_all[pav_all["gene_hs"].isin(genes_list)].rename(columns={"gene_hs": "gene"})
tax = pd.read_csv("tables/homology.species2group.csv")[["str", "group"]].rename(columns={"str": "species"})

# Filter groups
group_size = tax["group"].value_counts().reset_index()
group_size = group_size.rename(columns={"count": "species"})
group_size
group_select = group_size[group_size["species"] >= 10]["group"].tolist()
species_select = tax[tax["group"].isin(group_select)]["species"].tolist()
pav = pav[["gene"] + species_select]
pav.shape
group_size[group_size['group'].isin(group_select)].to_csv("tables/pav.group.selected.csv", index=False)
pav.to_csv("tables/pav.selected.csv", index=False)

# Deduplication
#pav_dedup = pav.groupby(species_select).apply(lambda x: x).reset_index(drop=True)
pav_dedup = pav.groupby(species_select)#.apply(lambda x: x).reset_index(drop=True)
pav_dedupped = pd.DataFrame()
pav_dedupped["gene_list"] = pav_dedup.apply(lambda x: ",".join(x["gene"]))
pav_dedupped["gene_num"] = pav_dedup.apply(lambda x: len(x))
pav_dedupped["gene_id"] = pav_dedup.apply(lambda x: ",".join(set(target).intersection(x["gene"])) if set(target).intersection(x["gene"]) else x["gene"].iloc[0])
pav_dedupped.reset_index(inplace=True)
pav_dedupped.shape
pav_deduped = pav_dedupped

pav_dd_mat = pav_deduped.set_index("gene_id")[species_select]#.to_numpy()
pav_dd_names = pav_deduped[["gene_id", "gene_num","gene_list",  ]]
#pav_dd_mat = pd.DataFrame(pav_dd_mat, index=pav_dd_names["gene_id"])
genes_select = pav_deduped["gene_id"].to_list()
pav_dd_names.to_csv("tables/pav.dedup.names.csv", index=False)


pav_dd_names = pd.read_csv("tables/pav.dedup.names.csv")

#%% Loss fraction
pav_dd_l = pav_deduped.melt(id_vars="gene_id", var_name="species", value_name="pav").merge(tax, on="species", how="left")
loss = pav_dd_l.groupby(["gene_id", "group"]).agg(loss=("pav", lambda x: sum(x == 0)), count=("pav", "count")).reset_index()
loss["frac"] = loss["loss"] / loss["count"]
loss_frac = loss.pivot(index="gene_id", columns="group", values="frac").reset_index()
lf_mat = loss_frac.set_index("gene_id")#.to_numpy()
lf_mat = lf_mat[group_select] # reorder columns
lf_mat.to_csv("tables/pav.dedupped.loss_frac.csv",index=True)

#%% Distance calculation
pav_g = pd.DataFrame({"group": group_select})
pav_g["mat"] = pav_g["group"].apply(lambda x: pav_dd_mat.loc[:, tax[tax["group"] == x]["species"]].to_numpy())
pav_g["d"] = pav_g["mat"].apply(lambda x: pdist(x, metric="cityblock") / x.shape[1])
pav_d_mat = np.vstack(pav_g["d"]).transpose()
pav_d_L2 = np.linalg.norm(pav_d_mat, axis=1)
pav_d = squareform(pav_d_L2)
# pav_gl = pav_g.explode("d").reset_index(drop=True)
# pav_gw = pav_gl.pivot(index="gene1", columns="group", values="dist").reset_index()
# pav_d = pav_gw.assign(d=lambda x: np.sqrt((x.select_dtypes(include=[np.number]) ** 2).sum(axis=1)))
pav_d.to_csv("tables/pav.dist.long.csv.gz", index=False)
pav_d_w = pav_d.pivot(index="gene1", columns="gene2", values="d").reset_index()
pav_d_w.to_csv("tables/pav.dist.wide.csv", index=False)
pav_d_mat = pav_d_w.set_index("gene1").to_numpy()
pav_d_df = pd.DataFrame(pav_d, index=genes_select, columns=genes_select)
pav_d_df.to_csv("tables/pav.dist.wide.py.csv", index=True)

#%% markers
# Assuming group_select and class_select are defined
group_select_df = pd.read_csv("tables/pav.group.selected.csv")
group_select = group_select_df['group'].tolist()
lf_mat = pd.read_csv("tables/pav.dedupped.loss_frac.csv")
lf_mat = lf_mat.set_index('gene_id')

# Create markers.1
markers_1 = np.eye(len(group_select))
markers_1_df = pd.DataFrame(markers_1, columns=group_select, index=group_select).reset_index().rename(columns={"index": "id"})
#print(markers_1_df)

# Create markers.2
markers_2 = pd.DataFrame([
    {"id": "NONE", "Aves": 0, "Fish": 0, "nrEutheria": 0, "Reptilia": 0, "Rodentia": 0},
    {"id": "Aves&Reptilia", "Aves": 1, "Fish": 0, "nrEutheria": 0, "Reptilia": 1, "Rodentia": 0},
    {"id": "Fish&Aves", "Aves": 1, "Fish": 1, "nrEutheria": 0, "Reptilia": 0, "Rodentia": 0},
    {"id": "Fish&Aves&Reptilia", "Aves": 1, "Fish": 1, "nrEutheria": 0, "Reptilia": 1, "Rodentia": 0},
    {"id": "ALL", "Aves": 1, "Fish": 1, "nrEutheria": 1, "Reptilia": 1, "Rodentia": 1}
])
markers_df = pd.concat([markers_1_df, markers_2], ignore_index=True)
markers_df = markers_df[['id'] + group_select]
print(markers_df)
markers_df.to_csv("tables/pav.marker.loss_frac.csv", index=False)
# Calculate distances
#lf_mat = lf.set_index("gene").to_numpy()
markers_df_d = markers_df.groupby("id").apply(lambda x: pd.DataFrame({
    "gene_id": lf_mat.index,
    "dist": np.linalg.norm(lf_mat - x.iloc[:, 1:].to_numpy(), axis=1)
    #"dist": sqrt(np.sum((lf_mat - x.iloc[:, 1:].to_numpy().reshape(1,-1))**2, axis=1))
}))
markers_df_d = markers_df_d.reset_index("id").reset_index(drop=True)


# Rank distances
markers_df_d["rank"] = markers_df_d.groupby("id")["dist"].rank(method="min")
markers_df_nn = markers_df_d[markers_df_d["rank"] <= 1]
print(markers_df_nn)


#%% rank markers
# Assuming markers_df_nn, target, pav_d_w, and markers_df_d are defined
markers_all = list(markers_df_nn["gene_id"]) + target
markers_labels = list(markers_df_nn["id"]) + target


tm_d = pav_d_df[markers_all].reset_index(names="gene_id")
tm_d_l = tm_d.melt(id_vars="gene_id", var_name="target", value_name="dist")
tm_d_l["rank"] = tm_d_l.groupby("target")["dist"].rank(method="min")


# set marker name
tm_d_l = tm_d_l.merge(markers_df_nn[["gene_id", "id"]].rename(columns={"gene_id":"target"}), on="target", how="left")
tm_d_l["marker"] = np.where(tm_d_l["id"].notna(), tm_d_l["id"], tm_d_l["target"])
tm_d_l = tm_d_l.drop(columns=["target", "id"])
tm_d_l

tm_d_l[tm_d_l['gene_id']==TM].sort_values("rank")


marker_neighbors = 30

markers_nb = markers_df_d[markers_df_d["rank"] <= marker_neighbors]

target_nb = tm_d_l[tm_d_l['marker'].isin(target)]
#target_nb = target_nb[target_nb['rank'] <= marker_neighbors + 1]
target_nb = target_nb[target_nb['rank'] <= marker_neighbors]
markers_nb = pd.concat([markers_nb.rename(columns={"id":"marker"}), target_nb])
markers_nb = markers_nb.merge(pav_dd_names, on="gene_id")
markers_nb = markers_nb.sort_values(["marker", "rank"])
markers_nb_col_reorder = ["marker", *[x for x in markers_nb.columns if x !="marker"]]
markers_nb = markers_nb[markers_nb_col_reorder]
markers_nb
markers_nb.to_csv("tables/pav.hmd.makers.v2.nb30.csv", index=False)


