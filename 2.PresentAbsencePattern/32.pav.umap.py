# generate umap plot coordinates
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
import umap
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA
from functools import partial
import itertools
import os

#%% clean data
# Read data
target = ["PRKAB2", "PRDM16", "DIO2", "UCP3", "PPARGC1A", "CIDEA", "UCP1", "ATP2A1", "ACTN3"]
TM = "TMEM41B"

group_select_df = pd.read_csv("tables/pav.group.selected.csv")
group_select = group_select_df["group"].tolist()


# # Gene loss type
lf_mat = pd.read_csv("tables/pav.dedupped.loss_frac.csv", index_col=0)
loss_frac_cut = 1/3
lf_l = lf_mat.reset_index().melt(id_vars="gene_id", var_name="group", value_name="loss_frac")
lf_l["loss_group"] = lf_l.apply(lambda x: x["group"] if x["loss_frac"] >= loss_frac_cut else "", axis=1)
gene_str = lf_l.groupby("gene_id").agg(loss_type=("loss_group", lambda x: "|".join(filter(None, x)))).reset_index()
gene_str["loss_type"] = gene_str["loss_type"].str.replace(r"\|+", "|").str.replace(r"^\||\|$", "")
gene_str_d = gene_str["loss_type"].value_counts().reset_index().rename(columns={"index": "loss_type", "loss_type": "num"})
gene_str["loss_str"] = gene_str["loss_type"].replace("", "NONE")

loss_type_dist = gene_str.loss_str.value_counts()
loss_type_dist
loss_type_most_freq = loss_type_dist.nlargest(7).index
gene_str["loss_name"] = gene_str['loss_str'].where(gene_str['loss_str'].isin(loss_type_most_freq), "Other")
gene_str['loss_name'] = gene_str['loss_name'].astype('category')


pav_d_df = pd.read_csv("tables/pav.dist.wide.py.csv", index_col=0)
pav_d_df


# #%% markers
# Create markers.1
markers_1 = np.eye(len(group_select))
markers_1_df = pd.DataFrame(markers_1, columns=group_select, index=group_select).reset_index().rename(columns={"index": "id"})
markers_2 = pd.DataFrame([
    {"id": "NONE", "Aves": 0, "Fish": 0, "nrEutheria": 0, "Reptilia": 0, "Rodentia": 0},
    {"id": "Aves&Reptilia", "Aves": 1, "Fish": 0, "nrEutheria": 0, "Reptilia": 1, "Rodentia": 0},
    #{"id": "Fish&Aves&Reptilia", "Aves": 1, "Fish": 1, "nrEutheria": 0, "Reptilia": 1, "Rodentia": 0},
    {"id": "ALL", "Aves": 1, "Fish": 1, "nrEutheria": 1, "Reptilia": 1, "Rodentia": 1}
])
markers_df = pd.concat([markers_1_df, markers_2], ignore_index=True)
markers_df = markers_df[['id']+group_select]

# Calculate distances
#lf_mat = lf.set_index("gene").to_numpy()
markers_df_d = markers_df.groupby("id").apply(lambda x: pd.DataFrame({
    "gene_id": lf_mat.index,
    "dist": np.sqrt(np.sum((lf_mat - x.iloc[:, 1:].to_numpy().reshape(1,-1))**2, axis=1))
    #'dist': (lf_mat - x.iloc[:, 1:].to_numpy()).to_list()
}))
markers_df_d = markers_df_d.reset_index("id").reset_index(drop=True)
markers_df_d.merge(lf_mat, left_on="gene_id", right_index=True, how="left")

# Rank distances
markers_df_d["rank"] = markers_df_d.groupby("id")["dist"].rank(method="min")
markers_df_nn = markers_df_d[markers_df_d["rank"] <= 1]
print(markers_df_nn)



#%% UMAP plot function
from plotnine import *
my_UMAP = partial(umap.umap_.UMAP, n_components=2,  metric="precomputed", 
                  n_jobs=1, random_state=2025 
                  )


#palette = ['#8DD3C7', '#FFFFB3', '#FB8072', '#80B1D3', '#FDB462', '#B3B3B3', '#BEBADA', '#FBBC3E'] #set3
palette = ['#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D', '#666666'] #dark2
palette_dict = dict(zip(gene_str['loss_name'].cat.categories, palette))

def umap_plot(n_nb, min_dist, spread, init_method, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    umap_config = my_UMAP( n_neighbors=n_nb, min_dist=min_dist, spread= spread, init=init_method)
    umap_result = umap_config.fit_transform(pav_d_df)
    base_name="{init}_s{spread:.1f}nb{n_nb}d{md:.2f}".format(
        init=umap_config.init, spread =umap_config.spread, n_nb=umap_config.n_neighbors, md=umap_config.min_dist )
    #base_name="d%.4fnb%ds%.1f" % (umap_config.min_dist, umap_config.n_neighbors, umap_config.spread)
    prefix = os.path.join(out_dir, base_name)


    umap_df = pd.DataFrame(umap_result, index=pav_d_df.index, columns=["V1", "V2"])
    umap_df.to_csv(prefix+".csv", index=True)
    #umap_df["focus"] = umap_df.index.isin(target)
    #umap_df["label"] = umap_df.apply(lambda x: x.name if x["focus"] else None, axis=1)
    umap_df = umap_df.merge(gene_str, left_index=True, right_on="gene_id", how="left")
    # markers_df_nn
    umap_target = umap_df[umap_df['gene_id'].isin([TM, *target])]
    umap_marker = umap_df[umap_df['gene_id'].isin(markers_df_nn['gene_id'])]
    umap_marker = umap_marker.merge(markers_df_nn[['id', 'gene_id']].rename(columns={"id":"marker"}))
    umap_target.loc[:,'marker'] = umap_target['gene_id']
    umap_highlight = pd.concat([umap_target, umap_marker])
    um_hl = umap_highlight
    

    p = (
        ggplot(umap_df, aes(x='V1', y='V2')) +
        geom_point(aes(color='loss_name'), shape='.') +
        #geom_path(aes(group='axis'), color="steelblue", data=um_axis) +
        geom_point(shape="o", data=um_hl,  size=1, color="black", alpha=0.5) +
        #geom_text_repel(data=um_hl, color="black") +
        geom_text(aes(label="marker"), size=8, data=um_hl, color="black") +
        scale_color_manual(values=palette_dict) +
        theme_bw() +
        theme(legend_position='bottom', legend_box_spacing=0) +
        coord_fixed() +
        labs(title=base_name )
    )
    # Save the plot
    p.save(prefix+".png", width=10, height=8)
    
    
# #%%% umap para search
# nb_vars = [15, 45]
# mdr_vars = [ 0.1, 1]
# spread_vars = [1, 5]
# init_vars = ["spectral"]
# 
# umap_para = pd.DataFrame(list(itertools.product(init_vars, spread_vars, nb_vars, mdr_vars)),
#                          columns=['init', 'spread', 'n_nb', 'min_dist_ratio'])
# #min_dist must be less than or equal to spread
# #umap_para = umap_para[umap_para['min_dist'] <= umap_para['spread']]
# umap_para['min_dist'] = umap_para['spread'] * umap_para['min_dist_ratio']
# for p in umap_para.itertuples(index=False):
#     print(p)
#     umap_plot(p.n_nb, p.min_dist, p.spread, p.init, out_dir="figures/hmd.umap")

#%% final plot    
nb = 15
min_dist = 1
spread = 1
init = "spectral"
umap_plot(nb, min_dist, spread, init, out_dir="figures/hmd.umap")







