# %%
import pandas as pd
import sys

file_path = sys.argv[1]

chromosomes = range(1, 27)


sheep_genome = { 1:278617202,
                2:250202058,
                3:226089100,
                4:121578099,
                5:108220788,
                6:118469697,
                7:101274418,
                8:91792871,
                9:95179658,
                10:86459471,
                11:62547497,
                12:80403655,
                13:83511835,
                14:66516657,
                15:82538637,
                16:71897364,
                17:73167223,
                18:67984460,
                19:60561550,
                20:51451717,
                21:47514194,
                22:51512491,
                23:62440644,
                24:42630236,
                25:44863754,
                26:45052359 }

sheep_genome_percentage = {chr: (length / sum(sheep_genome.values())) for chr,length in sheep_genome.items()}


for number, chromosome in enumerate(chromosomes):
    tmp_tajima = pd.read_csv(f'{file_path}/Tajima_{str(chromosome)}.csv')
    tmp_diversity = pd.read_csv(f'{file_path}/Diversity_{str(chromosome)}.csv')
    tmp_divergence = pd.read_csv(f'{file_path}/Divergence_{str(chromosome)}.csv')
    tmp_fst = pd.read_csv(f'{file_path}/Fst_{str(chromosome)}.csv')
    tmp_gnn = pd.read_csv(f'{file_path}/Gnn_{str(chromosome)}.csv', index_col = 0)
    if number == 0:
        combined_tajima = tmp_tajima
        combined_tajima.TajimaD = combined_tajima.TajimaD * sheep_genome_percentage[chromosome]
        combined_diversity = tmp_diversity
        combined_diversity.Diversity = combined_diversity.Diversity * sheep_genome_percentage[chromosome]
        combined_divergence = tmp_divergence
        combined_divergence = combined_divergence * sheep_genome_percentage[chromosome]
        combined_fst = tmp_fst
        combined_fst = combined_fst * sheep_genome_percentage[chromosome]
        combined_gnn = tmp_gnn
        combined_gnn = combined_gnn * sheep_genome_percentage[chromosome]
              
    else:
        combined_tajima.TajimaD = combined_tajima.TajimaD + (tmp_tajima.TajimaD * sheep_genome_percentage[chromosome])
        combined_diversity.Diversity = combined_diversity.Diversity + (tmp_diversity.Diversity * sheep_genome_percentage[chromosome])
        combined_divergence = combined_divergence + (tmp_divergence * sheep_genome_percentage[chromosome])
        combined_fst = combined_fst + (tmp_fst * sheep_genome_percentage[chromosome])
        combined_gnn = combined_gnn + (tmp_gnn * sheep_genome_percentage[chromosome])

combined_gnn.loc[:, "Breed"] = list(combined_gnn.index)

combined_tajima.to_csv("Combined_tajima.csv")
combined_diversity.to_csv("Combined_diversity.csv")
combined_divergence.to_csv("Combined_divergence.csv")
combined_gnn.to_csv("Combined_gnn.csv")
combined_fst.to_csv("Combined_fst.csv")

