# combine Ensembl results and convert copy number to gain/loss
library(tidyverse)
root_dir="data/homology/"
file_list = "data/homology/data.gene.list"
genes = scan(file_list, what=character())
genes
data_dir="data/homology/yuru_data/"
homo = tibble(gene=genes) %>% mutate(
  data = map(genes, ~read_csv(paste0(data_dir, .x, "_homologies.csv"), show_col_types = FALSE))
) %>% unnest(data)

dim(homo)
head(homo)
generic.skeleton()?scan
write_rds(homo, "data/homology/homology.RDS")
write_csv(homo, "data/homology/homology.all.csv.gz")

colnames(homo) = c("gene.hs", "homo_id", "method", "protein_id", "species", "tax_level", "homo_type")
species = homo %>% select(species) %>% distinct()
#species = species %>% mutate(Species = tools::toTitleCase(Species))
#species = species %>% mutate(Species = Species %>% str_replace("_", " "))
species
human = tibble(species="homo_sapiens")
species.all = rbind(species, human)

write_excel_csv(species.all, "data/homology/homo.species.csv")

species
homo.mm = homo %>% filter(Species =="mus_musculus")
write_excel_csv(homo.mm, "data/homology/homo.hs2mm.csv")


table(homo$method)
homo.pav = homo %>% select(gene.hs, species) %>% distinct()
homo.pav
pav.w = homo.pav %>% mutate(value=1) %>% 
  pivot_wider(names_from=species, values_from=value, values_fill=0)
pav.hs = tibble(gene.hs=pav.w$gene.hs, "homo_sapiens"=1)
pav.w.all = pav.w %>% left_join(pav.hs, by="gene.hs")
dim(pav.w.all)
write_excel_csv(pav.w.all, "data/homology/homology.pav.csv")
