library(tidyverse)
library(rtracklayer)
block_size = 45 
na_value = "-1"
TM="TMEM41B"

# prepare data ####
genes.list = read_csv("tables/pav.hmd.makers.v2.nb30.csv")
scores = read_delim("data/ensembl_canonical.gene.cds.phastCons.list.tsv", col_names=F)
colnames(scores) = c("gene_id", "length", "score_list")

## gene id to name 
## SLOW, the file should be prepared before the analysis
# prefix="data/Homo_sapiens.GRCh38.113"
# gtf_file = paste0(prefix, ".gtf.gz")
# gtf = import(gtf_file)
# gtf.df = gtf %>% as_tibble()
# ids = gtf.df %>% dplyr::select(type, tag, gene_name, gene_biotype, gene_id, transcript_id, protein_id)
# ids.q = ids %>% filter(!is.na(protein_id)) %>% distinct()
# write_excel_csv(ids.q, paste0(prefix, ".ids.csv"))
gene.id2name = read_csv("data/Homo_sapiens.GRCh38.113.ids.csv")
gene.id2name.raw = gene.id2name

scores = scores %>% separate(gene_id, into=c(NA, "id", NA), sep="\"", remove=F)
scores = scores %>% select(-gene_id) %>% rename(gene_id=id)
scores
genes.list.l = genes.list %>% mutate(gene = map(gene_list, ~unlist(strsplit(.x, ","))))
genes.list.l = genes.list.l %>% unnest(gene) 
marker.list = genes.list %>% pull(marker) %>% unique
dim(genes.list.l)
dim(genes.list)

gene.id2name = gene.id2name.raw %>% filter(tag=="Ensembl_canonical") %>% select(gene_name, gene_id) %>% distinct()
genes = genes.list.l %>% select(marker, gene)
scores.fmt = scores %>% left_join(gene.id2name, by=c("gene_id"))
scores.fmt
genes = genes %>% group_by(gene) %>% nest() %>% left_join(scores.fmt, by=c("gene"="gene_name"))
genes
# write_excel_csv(genes %>% unnest(data), "tables/phastCons.nb30.v2.csv")


# score to block ####
score2list = function(score_list) {
  s.df = tibble(value.str = strsplit(score_list, ",") %>% unlist())
  s.df = s.df %>% mutate(index=1:nrow(s.df), value=ifelse(value.str==na_value, NA, as.numeric(value.str)))
  s.df = s.df %>% mutate(block = ceiling(index/block_size))
  return(s.df)
}

score2block = function(score_df, block_size){
  score_df %>% mutate(block = ceiling(index/block_size)) %>%
   group_by(gene, marker, block) %>% summarise(
    value.sum = sum(value, na.rm=T),
    size = n(),
    na_num = sum(is.na(value)),
    nt = size - na_num,
    value.mean = ifelse(nt==0, NA, value.sum / nt),
    value.median = median(value, na.rm=T)
  )
}

# genes= read_csv("tables/phastCons.nb30.v2.csv")
marker.list = genes %>% pull(marker) %>% unique
genes = genes %>% mutate(score_block = map(score_list, score2list))
gene.score = genes %>%  unnest(score_block)
blocks = gene.score %>% score2block(block_size=block_size)
blocks = blocks %>% group_by(gene, marker) %>% mutate(block.frac=block/max(block))

# plot ####
plot_marker = function(marker_name, block_size=45, plot_genes=10, xlim_min=0) {
  score_df = gene.score %>% filter(marker == !!marker_name)
  blocks.t = score_df %>% score2block(block_size=block_size)
  blocks.t = blocks.t %>% group_by(gene, marker) %>% mutate(block.frac=block/max(block))
  blocks.tp = blocks.t %>% filter(nt >= block_size/2, !is.na(value.mean)) 
  blocks.tp
  blocks.tp.s = blocks.tp %>% group_by(gene, marker) %>% summarise(
    median=median(value.mean),
    Q1 = quantile(value.mean, 0.25),
    Q3 = quantile(value.mean, 0.75),
    IQR = Q3 - Q1,
    notch.low = Q1 - 1.5 * IQR,
    notch.high = Q3 + 1.5 * IQR)
  blocks.tp.s = blocks.tp.s %>%  arrange(desc(median)) %>% 
    ungroup %>% group_by(marker) %>% mutate(rank=rank(-median, ties.method="min"))
  blocks.tp.s
  blocks.tp.s.s = blocks.tp.s %>% group_by(marker) %>% summarise(marker.median=median(median))

  blocks.tp = blocks.tp %>% left_join(blocks.tp.s %>% select(gene, marker, median, rank, notch.low, notch.high))
  blocks.tp = blocks.tp %>% mutate(is.outlier = value.mean < notch.low | value.mean > notch.high)
  #blocks.tp = blocks.tp %>% arrange(rank) %>% mutate(gene.f = factor(gene, !!blocks.tp.s%>%pull(gene) %>% unique))
  blocks.tp = blocks.tp %>% arrange(rank) %>% mutate(gene.f = factor(gene, !!blocks.tp.s%>%pull(gene) %>% unique) )
  blocks.tp = blocks.tp %>% mutate(is.marker=ifelse(gene==!!TM,"TM", ifelse(gene %in% marker.list, "marker","NONE")))
  levels(blocks.tp$gene.f)
  blocks.tp %>% select(gene, marker) %>% distinct()

  blocks.tp.p = blocks.tp %>% filter((rank <= !!plot_genes)|(gene==marker))
  blocks.tp.p = blocks.tp.p %>% 
    mutate(rank.f=factor(rank, rev(seq(1, nrow(blocks.tp.p))))) %>%
    arrange(rank.f)
  blocks.tp.p = blocks.tp.p %>% mutate(value.mean=ifelse(value.mean < !!xlim_min, -Inf, value.mean))
  blocks.tp.p.names = with(blocks.tp.p, setNames(gene, rank.f))
  p = ggplot(blocks.tp.p , aes(x=value.mean, y=rank.f)) + facet_wrap(.~marker, scales="free_y") +
    geom_vline(aes(xintercept=v), color="lightgray", size=1/4, data=tibble(v=seq(0,1,0.1))) +
    geom_vline(aes(xintercept=marker.median), data=blocks.tp.s.s, color="dimgray", linetype="longdash") + #"darkviolet"
    geom_boxplot(aes(fill=is.marker), width=0.6, outliers=F, notch = TRUE) +
    geom_point(aes(color=block.frac), shape=1,  data=blocks.tp.p %>% filter(is.outlier), size=1.5) +
    scale_y_discrete(labels=blocks.tp.p.names) +
    scale_x_continuous(limits=c(xlim_min,1), breaks=seq(0,1, 0.2), expand=expansion(mult=0.01)) +
    scale_fill_manual("TMEM41B", values=c("TM"="green3", "marker"="gold", "NONE"="white")) +
    scale_color_gradient2("position", low="red", mid="darkgray", high="blue", midpoint = 0.5, 
                          limits=c(0, 1), breaks=c(0, 0.5, 1), labels=c("N-term", "mid", "C-term")) +
    labs(x="PhastCons", y="gene") +
    theme_classic() + theme(axis.title.y = element_blank()) + #legend.position = "none", 
    theme(panel.border = element_rect(fill="transparent", linewidth=1)) +
    guides(fill="none") 
  return(p)
}


## main figure ####
top.marker = "PPARGC1A"
p = plot_marker(top.marker, block_size=45, plot_genes = 6, xlim_min=0.4)
style = theme(panel.border = element_rect(fill="transparent", linewidth=1), axis.title.x = element_blank(),
              plot.margin=margin(r=5, unit="pt"))
p + style
ggsave("figures/conserve.phastCons.main.legend.svg", width=2.5, height=2, device="svg")



## plot all ####
marker.genes = c("PRKAB2", "PRDM16", "DIO2", "UCP3", "PPARGC1A", "CIDEA", "UCP1", "ATP2A1", "ACTN3")
marker.df = genes.list %>% select(marker) %>% distinct()
marker.df = marker.df %>% mutate(is.gene=marker %in% !! marker.genes)
marker.df = marker.df %>% arrange(desc(is.gene))
marker.df
# plots = marker.df %>% mutate(plot=map(marker, function(x)
#   plot_marker(x, block_size=45, plot_genes=6, xlim_min=0.4) + 
#     theme(axis.title.y = element_blank())))
# ggpubr::ggarrange(plotlist=plots$plot, common.legend=T, legend="bottom", align="v", ncol=4, nrow=5)
# ggsave("figures/conserve.phastCons.all.small.png", width=8, height=12)
# ggsave("figures/conserve.phastCons.all.small.pdf", width=8, height=12)

## Figure S1 ####
markers.psudo = marker.df %>% filter(!is.gene)
markers.psudo
plots.psudo = markers.psudo %>% mutate(plot=map(marker, function(x)
  plot_marker(x, block_size=45, plot_genes=5, xlim_min=0.4) + 
    theme(axis.title.y = element_blank())))
ggpubr::ggarrange(plotlist=plots.psudo$plot, common.legend=T, legend="bottom", align="v", ncol=3, nrow=4)
ggsave("figures/conserve.phastCons.psudo.small.png", width=8, height=8)
ggsave("figures/conserve.phastCons.psudo.small.pdf", width=8, height=8)


