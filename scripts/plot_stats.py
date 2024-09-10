# %%
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from plotnine import *
import seaborn as sns
import colorcet as cc 
import pandas as pd
import numpy as np
import scipy
import os


# %% [markdown]
# Time to most recent common ancestor

# %%
fig, ax = plt.subplots(figsize=(15, 7))    
ax.step(tmrca_df['window'], tmrca_df['time'], color='black', linewidth=0.5)
ax.set_ylabel('TMRCA (generations)')

fig.savefig('TMRCA.png')

# %% [markdown]
# GNN

# %%
gnn_df = pd.read_csv(f"Combined_gnn.csv", index_col=0)
# breeds=list(set(gnn_df.Breed))

last_col_index = list(gnn_df.columns).index('Breed')
# Keep only the columns from 0 to "population"
gnn_df = gnn_df.iloc[:, :last_col_index + 1]
gnn_df = gnn_df.groupby('Breed').mean()
# zscore normalise
for col in list(gnn_df):
    gnn_df[col] = np.log(gnn_df[col])
    # gnn_df[col] = scipy.stats.log(gnn_df[col])
    
row_linkage = scipy.cluster.hierarchy.linkage(gnn_df, method='average')
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = gnn_df.index.values[order]

palette_s = sns.color_palette(cc.glasbey, n_colors=len(x_pop))
subsp_clr = {x_pop[i]: clr for i, clr in enumerate(palette_s)}

colours = pd.Series(subsp_clr) 
cg = sns.clustermap(gnn_df[x_pop], row_linkage=row_linkage, col_cluster=False,
                    figsize=(10, 10), rasterized=True, 
                    #method='ward', robust=True, 
                    row_colors=colours)
cg.ax_heatmap.set_yticks([])
cg.ax_heatmap.set_xticks([])
# cg.ax_heatmap.set_ylabel(label)
cg.cax.remove()

# Adjust these values to move the colorbar
cbar_ax = cg.fig.add_axes([.1, .07, .15, .015])  
plt.colorbar(cg.ax_heatmap.get_children()[0], cax=cbar_ax, orientation='horizontal')
for pop, col in colours.items():
    cg.ax_col_dendrogram.bar(0, 0, color=col, label=pop, linewidth=0)

# Create custom legend handles and labels
handles = [mpatches.Patch(color=color, label=label) for label, color in subsp_clr.items()] 
cg.ax_col_dendrogram.legend(handles=handles, 
                            ncol=len(subsp_clr)/2, 
                            prop={'size': 12}, 
                            loc='upper center', bbox_to_anchor=(0.5, 0.8)
                            )      
cg.figure.savefig(f"GNN.png", bbox_inches='tight', transparent=True, dpi=300)
plt.close()

# %% [markdown]
# Non-windowed stats

# %%
tajima = pd.read_csv("Combined_tajima.csv")
tajima = tajima.sort_values(by='TajimaD').reset_index(drop = True)
tajima['Breed'] = pd.Categorical(tajima.Breed, categories=pd.unique(tajima.Breed))

# %%
p = ggplot(data = tajima) + \
    geom_bar(aes(x = "reorder(Breed, TajimaD)", y = "TajimaD", fill='Breed'), stat = "identity") + \
        ggtitle(f"TajimaD, not windowed")+ \
            scale_fill_manual(values=["#4D4887", "#77E6B0"] * len(set((tajima.Breed)))) + \
                theme(legend_position = "none", 
                    axis_text_x = element_text(angle = 90),
                    axis_title_x = element_blank())

# %%
ggsave(p, f"TajimaD.png")

# %%
diversity = pd.read_csv("Combined_diversity.csv")

# %%
diversity = diversity.sort_values(by='Diversity').reset_index(drop = True)
diversity['Breed'] = pd.Categorical(diversity.Breed, categories=pd.unique(diversity.Breed))

# %%
p = ggplot(data = diversity) + \
    geom_bar(aes(x = "Breed", y = "Diversity", fill='Breed'), stat = "identity") + \
        ggtitle(f"Diversity, not windowed") + \
            scale_fill_manual(values=["#4D4887", "#77E6B0"] * len(set((diversity.Breed)))) + \
                theme(legend_position = "none", 
                    axis_text_x = element_text(angle = 90),
                    axis_title_x = element_blank())

# %%
ggsave(p, f"Diversity.png")

# %%
Fst = pd.read_csv("Combined_fst.csv", index_col=0)

# %%
# add breeds as index
Fst.index = Fst.columns

# Mask the upper triangle and diagonal
mask = np.triu(np.ones_like(Fst, dtype=bool))
# Remove the first row and last column from the DataFrame
Fst = Fst.iloc[1:, :-1]
mask = mask[1:, :-1]
# Plot as a heatmap
fig, ax = plt.subplots(figsize=(6, 6))
cax = ax.matshow(Fst.mask(mask).astype(float), cmap='viridis')
# Add colorbar
fig.colorbar(cax)
# Adjust ticks and labels
ax.set_xticks(range(len(Fst.columns)))
ax.set_yticks(range(len(Fst.index)))
ax.set_xticklabels(Fst.columns, rotation=90)
ax.set_yticklabels(Fst.index)
# Adjust layout
plt.tight_layout()

fig.savefig('Fst_heatmap.png', dpi=300) 


divergence = pd.read_csv("Combined_divergence.csv", index_col=0)

# %%
# add breeds as index
divergence.index = divergence.columns

# Mask the upper triangle and diagonal
mask = np.triu(np.ones_like(divergence, dtype=bool))
# Remove the first row and last column from the DataFrame
divergence = divergence.iloc[1:, :-1]
mask = mask[1:, :-1]
# Plot as a heatmap
fig, ax = plt.subplots(figsize=(6, 6))
cax = ax.matshow(divergence.mask(mask).astype(float), cmap='viridis')
# Add colorbar
fig.colorbar(cax)
# Adjust ticks and labels
ax.set_xticks(range(len(divergence.columns)))
ax.set_yticks(range(len(divergence.index)))
ax.set_xticklabels(divergence.columns, rotation=90)
ax.set_yticklabels(divergence.index)
# Adjust layout
plt.tight_layout()

fig.savefig('Divergence_heatmap.png', dpi=300) 
