library(rtracklayer)
library(tidyverse)

gtf_file = "data/Homo_sapiens.GRCh38.113.gtf.gz"
gtf = import(gtf_file)
gtf.df = gtf %>% as_tibble()
colnames(gtf.df)

# table(gtf$type)
gtf.gene = gtf %>% as_tibble() %>% filter(type=="gene") %>%
  select(source, gene_id, gene_name, gene_biotype, tag)


gtf.gene
gtf.gene.stat = gtf.gene %>% group_by(source, gene_biotype, tag) %>% summarise(num=n())
gtf.gene.stat = gtf.gene.stat %>% ungroup() %>% arrange(source, desc(num))
gtf.gene.stat
# 
gtf.transcript = gtf %>% as_tibble() %>% filter(type=="transcript") %>%
  select(source, gene_id, gene_name, gene_biotype,
                                          tag, transcript_id, transcript_name, transcript_biotype)
gtf.transcript.summary = gtf.transcript %>%
  group_by(source, gene_biotype, transcript_biotype, tag) %>% summarise(num=n())
gtf.transcript.summary = gtf.transcript.summary %>% ungroup() %>% arrange(source, desc(num))
gtf.transcript.summary

gtf.transcript.coding = gtf.transcript %>% 
  filter(gene_biotype=="protein_coding") 
  # do not filter transcript_biotype, as the Ensembl_canonical transcript for protein coding gene 
  # may not be a protein_coding transcript, eg. for gene PKD1L2
  # transcript_biotype=="protein_coding"

coding = gtf.transcript.coding

coding.gene = coding %>% group_by(gene_id, tag) %>% summarise(num=n())
coding.gene.w = coding.gene %>% pivot_wider(id_cols = gene_id, names_from=tag, values_from=num, values_fill=0)
coding.gene.w
table(coding.gene.w$Ensembl_canonical)

# 
# coding.gene.w %>% filter(Ensembl_canonical==0)
# coding %>% filter(gene_id == "ENSG00000166473")
# gtf.df %>%  filter(gene_id == "ENSG00000166473")
# tmp = .Last.value

table(gtf.df$tag)
table(coding$tag)
table(gtf.df$type)


canonical = coding %>% filter(tag =="Ensembl_canonical") %>% pull(transcript_id)
coding.cds.canonical = subset(gtf, type=="CDS" & transcript_id %in% canonical)

coding.cds.canonical %>% as_tibble()
table(coding.cds.canonical%>% as_tibble() %>% pull(tag)) 

coding.cds = subset(gtf, gene_biotype=="protein_coding" & type=="CDS" & tag=="Ensembl_canonical")
table(coding.cds %>% as_tibble() %>% pull(tag)) 

coding.cds.canonical
result.clean = coding.cds.canonical
colnames(mcols(coding.cds.canonical))
mcols(result.clean) <- mcols(result.clean)[, c("source", "type", "score", "phase", "gene_id")]
result.clean
export(result.clean, con="data/ensembl_canonical.cds.clean.ensembl.gtf", format="GTF")



# mapping the ensembl chrom name to phastCons NCBI convention
chr.map = read_csv("data/hg38.chrom.convert_table.csv")
chr.map$chr.ensembl
seqlevels(gtf)
chr.mapping = setNames(chr.map$chr.ncbi, chr.map$chr.ensembl)
seqlevels(coding.cds.canonical) 
seqlevels(coding.cds.canonical) <- chr.mapping[seqlevels(coding.cds.canonical)]
coding.cds.canonical
export(coding.cds.canonical, con="data/ensembl_canonical.cds.gtf", format="GTF")
write_excel_csv(coding.cds.canonical %>% as_tibble(), "data/ensembl_canonical.cds.csv")

result.clean = coding.cds.canonical
colnames(mcols(coding.cds.canonical))
mcols(result.clean) <- mcols(result.clean)[, c("source", "type", "score", "phase", "gene_id")]
result.clean
export(result.clean, con="data/ensembl_canonical.cds.clean.ncbi.gtf", format="GTF")

