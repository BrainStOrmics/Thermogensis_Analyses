library(tidyverse)
species.df = read_csv("tables/homology.species2group.csv")
species.vec = species.df %>% pull(str)
#align = read_csv("tables/seq.pair_align.p5sim.csv")
#align = read_csv("tables/seq.pair_align.dm67.p5sim.csv")
align = read_csv("tables/seq.pair_align.hs67.p5sim.csv")
align.raw = align
align = align.raw %>% filter(source.species %in% !! species.vec, target.species %in% !! species.vec)
dim(align)
align = align %>% group_by(source.species, target.species) %>% 
  arrange(p3d, p5d) %>% mutate(rank=row_number())
align = align %>% ungroup() %>% filter(rank==1)
dim(align)
align.species = c(align %>% select(source.species) %>% distinct() %>% pull(source.species),
                  align %>% select(target.species) %>% distinct() %>% pull(target.species)
                  ) %>% unique()
length(align.species)
196*196 - 196
species.vec[!(species.vec %in% align.species)]

align.full = expand_grid(source.species=align.species, target.species=align.species)
align.full = align.full %>% left_join(align)
align.miss = align.full %>% filter(is.na(p5d), source.species!=target.species)
align.miss
# 1 pair align missed. ignore it

align.l = align %>% select(source.species, target.species, p5d, p3d) %>% pivot_longer(c(p5d, p3d))
align.l = align.l %>% mutate(name.f = factor(name, levels=c('p5d', 'p3d'), labels=c("N-term", "rest")))
ggplot(align.l, aes(x=name.f, y=value, fill=name.f)) +
  #geom_violin(alpha=0.5, color="gray") +
  #geom_boxplot(aes(color=name.f),notch = F, outliers = F, coef=0, fill="transparent", width=0.4) +
  geom_violin(alpha=1) + 
  geom_segment(x=1, xend=2, y=1.05, data=tibble(x=1), inherit.aes = F) +
  geom_text(x=1.5, y=1.05, label="****", data=tibble(x=1), inherit.aes = F, vjust=0) +
  theme_classic() + 
  labs(x="region", y="pairwise normalized Hamming distance") +
  scale_y_continuous(limits=c(0, 1.1)) +
  #scale_color_manual(values=c("red", "blue")) +
  scale_fill_manual(values=c("maroon", "darkgreen")) +
  theme(legend.position = "none")
#ggsave("figures/seq.pair_align.dm67.p5sim.png", width=3, height=4)
#ggsave("figures/seq.pair_align.dm67.p5sim.svg", width=3, height=4, device="svg")
ggsave("figures/seq.pair_align.hs67.p5sim.png", width=3, height=4)
ggsave("figures/seq.pair_align.hs67.p5sim.svg", width=3, height=4, device="svg")


wilcox.test(align$p5d, align$p3d, paired=T)
t.test(align$p5d, align$p3d, paired = T)


