library(tidyverse)
library(ggside)
library(Hmisc)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(ggrepel)
library(aplot)

# get top freq genes ####
# final version
df = read_delim("tables/gene_ranking_up36down31.txt")


#df = df %>% rename(up_count=Up, down_count=Down)
dfd = df %>% group_by(up_count, down_count) %>% summarise(num=n()) %>% ungroup()



up_q = df %>% filter(up_count > 0) %>% pull(up_count) %>% quantile(0.9) %>% floor
down_q = df %>% filter(down_count > 0) %>% pull(down_count) %>% quantile(0.9) %>% floor
table(df$up_count)
up_q
down_q


# write top genes
df.selected = df %>% filter(up_count >= up_q | down_count >= down_q)

# plot ####
up_limit =max(df$up_count) +1
down_limit = max(df$down_count) + 1
axis_limit = max(up_limit, down_limit)
axis_limit
down_limit = up_limit
frac_limit = 0.23
dfd
## major heatmap ####
p =ggplot(dfd, aes(x=up_count, y=down_count)) +
  #geom_point(size=1) +
  geom_tile( aes(fill=log2(num))) +
  geom_polygon(aes(x=x, y=y),fill="black", alpha=0.7, inherit.aes = F,
            data=tibble(x=c(-Inf, -Inf, up_q, up_q), y=c(-Inf, down_q, down_q, -Inf))) +
  geom_hline(yintercept=down_q, color="red", linetype="dashed") +
  geom_vline(xintercept=up_q, color="red", linetype="dashed") +
  scale_y_continuous(limits=c(-0.5, axis_limit), expand=expansion()) +
  scale_x_continuous(limits=c(-0.5,axis_limit), expand=expansion()) +
  labs(x="Up frequency", y="Down frequency") +
  scale_color_distiller(palette="Spectral", aesthetics = "fill", breaks=seq(0,10,2)) +
  #coord_fixed() +
  #theme_bw()
  theme_classic()

p
p = p + guides(fill = guide_colorbar(barwidth = 4, barheight = 1, title.position = "top", title.hjust = 0.5)) +
  theme(legend.position = "bottom")
p
legend.p = get_legend(p)


## up freq margin ####
df.up = df %>% filter(up_count > 0)
quantile(df.up$up_count, 0.9)
up.cumu = tibble(up_count=seq(0, axis_limit, 1)) %>% mutate(
  cumu=map(up_count, ~sum(df.up$up_count>=.x) /length(df.up$up_count))
) %>% unnest(cumu)

up.p <- ggplot(up.cumu, aes(x=up_count, y=cumu)) +
  geom_line(color="blue") +
  geom_ribbon(aes(xmin=-Inf, xmax=up_count, y=cumu), fill="blue", alpha=0.3) +
  geom_segment(x=-Inf, xend=up_q, y=0.1, color='red', linewidth=1/2, 
               inherit.aes = F, data=tibble(row=1)) +
  geom_segment(x=up_q, y=0.1, yend=-Inf, color="red", linewidth=1/2, linetype="dashed", inherit.aes = F,
               data = tibble(row=1)) +
  #theme_bw() +
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_y_continuous("%", minor_breaks = seq(0,1,0.1), #breaks=c(0,0.1, 0.5,1), 
                     breaks=c(0, 0.1, 0.2),
                     expand=expansion(),
                     labels=function(x)x*100) + #labels=scales::percent
  scale_x_continuous(limits=c(-0.5, axis_limit), expand=expansion(0)) +
  coord_cartesian(ylim=c(0, frac_limit))
  #theme(plot.background = element_blank())
up.p

## down freq margin ####
df.down = df %>% filter(down_count > 0)
down.cumu = tibble(down_count=seq(0, axis_limit, 1)) %>% mutate(
  cumu=map(down_count, ~sum(df.up$down_count>=.x) /length(df.up$down_count))
) %>% unnest(cumu)

down.p <- 
  ggplot(down.cumu, aes(y=down_count, x=cumu)) +
  geom_line(color="darkgreen") +
  geom_ribbon(aes(ymin=-Inf, ymax=down_count), fill="darkgreen", alpha=0.3) +
  geom_segment(y=-Inf, yend=down_q, x=0.1, color='red', linewidth=1/2, 
               inherit.aes = F, data =tibble(row=1)) +
  geom_segment(y=down_q,x=0.1, xend=-Inf, color="red", linewidth=1/2, linetype="dashed", 
               inherit.aes = F, data = tibble(row=1)) +
  #theme_bw() +
  theme_classic() +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  scale_x_continuous("%",  minor_breaks = seq(0,1,0.1),
                     breaks=c(0,0.1, 0.2, 0.5,1),
                     #limits=c(0, frac_limit),
                     expand=expansion(),
                     labels=function(x)x*100) +
  scale_y_continuous(expand=expansion(), limit=c(-0.5, axis_limit)) +
  coord_cartesian(xlim=c(0,frac_limit))
down.p

## legends ####
empty.p <- ggplot() + theme_void() 

# plot_grid(up.p, NULL, p + theme(legend.position = "none"), down.p,
#         ncol=2, align="hv", axis="lb", rel_heights = c(1,2), rel_widths = c(2,1))
legend.p <- ggplot(tibble(x=1)) +
  geom_rect(xmin=0, xmax=1, ymin=0.2, ymax=0.8, fill="white", color="white", alpha=0.7) +
  #geom_text(label="filtered", x=0.5, y=0.5,color="white", size=3) +
  geom_label(label="filtered", x=0.5, y=0.5, fill="black", alpha=0.7,color="white", size=3,
              label.r=unit(0,"pt")) +
  theme_void()
legend.p

## combine plots ####
df
known.plot = c("UCP1", "UCP3", "PPARGC1A")
df.marker = df %>% filter(gene %in% !!known.plot)
p.labeled <- p + geom_point(data=df.marker,shape=16, size=1) +
  geom_text_repel(aes(label=gene), data=df.marker)

p.combined <- (p.labeled + theme(legend.position = "none"))  %>% 
  aplot::insert_top(up.p, 1/4) %>%
  aplot::insert_right(down.p, 1/4)

p.combined
# aplot is not compatible with coord_fixed
#p.combined[2,1] = p.combined[2,1] + coord_fixed()

ggsave("figures/up_down_freq.svg", width=3, height=3)

