# %% [markdown]
# # Population statistics
# 
# try to do some statistics on smarter samples

# %%
import tskit
import json
import pandas as pd
import numpy as np
import sys

from collections import Counter
from tskitetude import get_project_dir

# %% [markdown]
# Determine chromosome and file to open

# %%
chromosome = sys.argv[1]
file_path = sys.argv[2]
ts_name = sys.argv[3]

tsFile = f"{file_path}/{ts_name}.trees"

# %%
ts = tskit.load(tsFile)
tsPos = [x.position for x in ts.sites()]

# %% [markdown]
# Preapre some lists

# %%
breeds = list(set([json.loads(ts.population(ts.node(u).population).metadata)['breed'] for u in ts.samples()]))

sample_nodes = [ts.node(n) for n in ts.samples()]
samples_listed_by_breed_dict = { pop: [s.id for s in sample_nodes if json.loads(ts.population(s.population).metadata)['breed'] == pop] for pop in breeds}

samples_listed_by_breed = [ [s.id for s in sample_nodes if json.loads(ts.population(s.population).metadata)['breed'] == pop] for pop in breeds]

# %%
num_populations = ts.num_populations
breedPairs = [(x, y) for x in range(num_populations) for y in range(num_populations) if x < y]

windows = np.linspace(0, int(ts.sequence_length), num= int(ts.sequence_length) // 100)

# %% [markdown]
# Time to most recent common ancestor

# %%
tmrca = ts.divergence(sample_sets=samples_listed_by_breed, indexes=breedPairs, windows=windows, mode='branch') // 2

# %% [markdown]
# Genealogical nearest neighbor

# %%
gnn = ts.genealogical_nearest_neighbours(
    ts.samples(), samples_listed_by_breed
)

cols = {breed: gnn[:, u] for u, breed in enumerate(breeds)}
# cols["breed"] = [json.loads(ts.population(ts.node(u).population).metadata)["breed"] for u in ts.samples()]
GnnDF = pd.DataFrame(cols, index=[json.loads(ts.population(ts.node(u).population).metadata)["breed"] for u in ts.samples()])
GnnDF.to_csv(f"Gnn_{chromosome}.csv")

# %% [markdown]
# Non-Windowed stats

# %% [markdown]
# Tajima by breed

# %%
popsTajima = pd.DataFrame()

for breed in breeds:
    tajima = ts.Tajimas_D(sample_sets=samples_listed_by_breed_dict[breed])
    tmp = pd.DataFrame({"TajimaD": [tajima], "Breed": breed})
    popsTajima = pd.concat([popsTajima, tmp])

popsTajima.to_csv("Tajima_" + str(chromosome) + ".csv", index = None)


# %% [markdown]
# Diversity by breed

# %%
popsDiversity = pd.DataFrame()
for breed in breeds:
    diversity = ts.diversity(sample_sets=samples_listed_by_breed_dict[breed])
    tmp = pd.DataFrame({"Diversity": [diversity], "Breed": breed})
    popsDiversity = pd.concat([popsDiversity, tmp])

popsDiversity.to_csv("Diversity_" +str(chromosome) + ".csv", index = None)

# %% [markdown]
# Fst by breed

# %%
Fst = ts.Fst(samples_listed_by_breed, indexes=breedPairs, windows=None, mode='site', span_normalise=False)

# fst = np.average(fst_values, axis=0, weights=chr_lengths)
fst_df = pd.DataFrame(index=breeds, columns=breeds)
for (i, j), fst_value in zip(breedPairs, Fst):
    fst_df.iloc[i, j] = fst_value
    fst_df.iloc[j, i] = fst_value

fst_df.to_csv("Fst_" + str(chromosome) + ".csv", index = None)    

# %%
# Divergence by breed pairs
divergence = ts.divergence(samples_listed_by_breed, indexes=breedPairs, windows=None, mode='site', span_normalise=False)

divergence_df = pd.DataFrame(index=breeds, columns=breeds)
for (i, j), values in zip(breedPairs, divergence):
    divergence_df.iloc[i, j] = values
    divergence_df.iloc[j, i] = values

divergence_df.to_csv("Divergence_" + str(chromosome) + ".csv", index = None) 

