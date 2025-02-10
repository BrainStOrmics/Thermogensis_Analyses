# refine umap plot
library(tidyverse)
library(umap)
library(Rtsne)
library(ggrepel)
library(stats)
library(vegan)
library(ggdendro)
library(colorspace)


target = c("PRKAB2", "PRDM16", "DIO2", "UCP3", "PPARGC1A", "CIDEA", "UCP1", "ATP2A1", "ACTN3")
datadir = "figures/hmd.umap/"
prefix="spectral_s1.0nb15d1.00"
# read data
#umap_file = paste0(datadir, prefix, ".csv")
umap_file = "tables/hmd.umap.spectral_s1.0nb15d1.00.csv"
um_df = read_csv(umap_file)
colnames(um_df) = c("gene_id", "UMAP1", "UMAP2")
dim(um_df)

loss.frac = read_csv("tables/pav.dedupped.loss_frac.csv")
groups.df = tibble(group=c("nrEutheria", "Fish", "Rodentia", "Reptilia", "Aves"))
#groups.df = read_csv("tables/pav.group.selected.csv")
#markers.df = read_csv("tables/pav.hmd.makers.nb30.csv")
markers.df = read_csv("tables/pav.hmd.makers.v2.nb30.csv")
groups.df

## gene loss.type ####
loss.frac.cut = 0.5
loss.frac.bin.break = c(0, 5,25,50,100)/100
lf.l = loss.frac %>% pivot_longer(-gene_id, names_to="group", values_to="loss.frac")
lf.l = lf.l %>% mutate(loss.group=ifelse(loss.frac >= loss.frac.cut, group, ""))
gene.str = lf.l %>% group_by(gene_id) %>% 
  summarise(loss.type=paste(loss.group, collapse="|"),
            loss.groups = sum(loss.group != ""),
            major.frac = max(loss.frac),
            major.group = group[which.max(loss.frac)])
gene.str = gene.str %>% mutate(loss.type=loss.type %>% str_replace_all("\\|+","|") %>% str_replace_all("^\\||\\|$", ""))
gene.str = gene.str %>% mutate(loss.str=ifelse(loss.type=="", "NONE", loss.type))
gene.str

gene.str = gene.str %>% mutate(major.frac.bin = cut(major.frac, loss.frac.bin.break, include.lowest = T, right = F))
gene.str = gene.str %>% mutate(pattern=ifelse(major.frac < loss.frac.bin.break[2], "NONE", ifelse(loss.groups<=1, major.group, loss.str)))
gene.str.d = gene.str %>% group_by(pattern, major.frac.bin) %>% summarise(num=n()) %>% arrange(desc(num))
gene.str.d
table(gene.str$loss.type) %>% enframe("loss.type", "num") %>% arrange(desc(num))
 
gene.patterns = c("NONE", groups.df$group,"Fish|Reptilia|Aves", "Reptilia|Aves", "Fish|Aves")

gene.str = gene.str %>% mutate(pattern = factor(pattern) %>% fct_infreq %>% fct_other(keep=gene.patterns))
gene.str
#gene.str = gene.str %>% mutate(loss.type=factor(loss.str) %>% fct_infreq %>% fct_lump(7))
table(gene.str$loss.str) %>% enframe() %>% arrange(desc(value))
#write_excel_csv(gene.str, "tables/pav.loss_pattern.csv")


# markers ####
markers.clean = markers.df %>% filter(rank<=1) %>% select(marker, gene_id, rank) 
markers.clean = bind_rows(markers.clean, tibble_row(marker="TMEM41B", gene_id="TMEM41B", rank=1))
um.markers = markers.clean %>% left_join(um_df)
um.hl = um.markers# %>% filter(marker != "Fish&Aves&Reptilia")
um.hl = um.hl %>% mutate(is.gene=marker %in% !!target)
um.hl %>% filter(!is.gene)
um.label = c("NONE"="Present\nin all", "ALL"="Absent\nin all",
             "Aves"="Absent in\nAves", "Fish"="Absent in\nFish",
             "Reptilia"="Absent in\nReptilia",
             "Rodentia"="Absent in\nRodentia",
             "nrEutheria"="Absent in\nnrEutheria",
             "Aves&Reptilia"="Absent in\nSauropsida",
             "Fish&Aves"="Absent in\nFish&Aves",
             "Fish&Aves&Reptilia"="Absent in\nFish&Sauropsida"
             )

um.hl = um.hl %>% left_join(um.label %>% enframe("marker", "label"), by="marker")
um.hl = um.hl %>% mutate(label = ifelse(is.na(label), marker, label))
um.hl



# plot config ####
single_color = setNames(scales::hue_pal()(5), c("nrEutheria", "Fish", "Reptilia", "Aves", "Rodentia"))
hybrid_color=c("Reptilia|Aves"="darkgreen", "Fish|Aves"="darkred", "Fish|Reptilia|Aves"="darkblue")
colors = c("NONE"="gray",single_color, hybrid_color,  "Other"="darkgray")
colors
colors.v2 = c("NONE"="darkgray", single_color, hybrid_color, "Other"="black")
colors.v3 = colorspace::desaturate(colors.v2, 0.5)
#RColorBrewer::display.brewer.all()
pattern.vec = names(colors)
pattern.vec

label.manual=c("nrEutheria",  # top left
               "Rodentia", "Aves", "Fish&Aves", "ALL", # left
               "PRDM16", "Aves&Reptilia", "ATP2A1", "Fish&Aves&Reptilia",  # manual
               "Fish", "UCP1", # bottom right
               "UCP3", "ACTN3" #top right
               )
um.hl.m = um.hl %>% filter(marker %in% label.manual)
text.size =4
nudge=0.2

# plot ####
um.df = um_df %>% left_join(gene.str%>% select(gene_id, loss.groups, pattern, major.frac.bin,), by="gene_id")
um.df = um.df %>% mutate(pattern=factor(pattern, pattern.vec))
ggplot(um.df, aes(x=UMAP1, y=UMAP2)) + 
  geom_point(aes(color=pattern, fill=pattern, shape=major.frac.bin, alpha=major.frac.bin), 
             data=um.df %>% filter(pattern!="NONE")) +
  geom_point(aes(color=pattern), shape=16, data=um.df %>% filter(pattern=="NONE")) +
  geom_point(data=um.hl, shape=3, size=2, color="black") +
  geom_text_repel(aes(label=label), data=um.hl, color="black", size=4,
                  lineheight = 0.8) +
  scale_color_manual( values = colors.v3,
                     aesthetics = c("fill", "color"),
                     breaks=c("NONE", "Fish", "Aves", "Reptilia", "Rodentia", "nrEutheria", 
                              names(hybrid_color), "Other"),
                     labels=function(x)ifelse(x=="NONE", "AllPresent",
                                              gsub("Reptilia\\|Aves","Sauropsida", x ))
  ) +
  scale_alpha_manual(values=c(0.4, 0.4, 1) %>% rev,labels=c("(0,0.25)","[0.25,0.5)", "[0.5,1]")) +
  scale_shape_manual(values=c(1, 16, 16) %>% rev, labels=c("(0,0.25)","[0.25,0.5)", "[0.5,1]")) +
  theme_bw() + 
  coord_fixed() + 
  labs(color="Absent\nPattern", fill="Absent\nPattern", 
       alpha="Absent\nProportion", shape="Absent\nProportion") #+

version = "v6"
file_prefix = paste0("figures/umap.", prefix, ".", version)
ggsave(paste0(file_prefix, ".png"), width=10, height=6)
ggsave(paste0(file_prefix, ".svg"), width=10, height=6, device="svg")


