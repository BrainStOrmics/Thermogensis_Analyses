library("tidyverse")
species = read_csv("data/homology/homo.species.csv")
species

genus.tax = read_delim("data/taxonomy/homo.genus.tax.txt", delim="\t")
genus.tax
colnames(genus.tax)
taxonomy.levels <- c(
  "Superkingdom",
  "Kingdom",
  "Subkingdom",
  "Phylum",
  "Subphylum",
  "Superclass",
  "Class",
  "Subclass",
  "Infraclass",
  "Cohort",
  "Subcohort",
  "Superorder",
  "Order",
  "Suborder",
  "Infraorder",
  "Parvorder",
  "Superfamily",
  "Family",
  "Subfamily",
  "Tribe",
  "Subtribe",
  "Genus",
  "Subgenus",
  "Species"
) %>% tolower()
genus.tax.ordered = genus.tax[,c("name", taxonomy.levels, "tax_str")]
genus.tax.l = genus.tax %>%  select(-tax_str) %>% pivot_longer(-name, names_to="level")
genus.tax.l

species = species %>% separate(species, c("name", NA, NA), sep="_", remove=F)
species
tax = species %>% rename(str=species) %>% left_join(genus.tax.ordered, by="name")
#write_excel_csv(tax, "data/homology/homo.species.tax.csv")

tax
tax.l = tax %>% select(-name, -tax_str) %>% pivot_longer(-str, names_to="level")
level.stat = table(tax.l$level) %>% enframe()
level.stat
level.stat$name[!(level.stat$name %in% taxonomy.levels)]
tax.l = tax.l %>% mutate(level = factor(level, levels=taxonomy.levels))

tax.l.s = tax.l %>% group_by(level) %>% summarise(uniq.num=length(unique(value)), na.num=sum(is.na(value)), num=n())

levels.base = c("phylum", "class", "order", "family")
tax.l.s %>% filter(level %in% levels.base)

tax.stat = tax %>% group_by(phylum, class, order, family) %>% summarise(num=n())
#write_excel_csv(tax.stat, "data/homology/homo.species.tax2group.csv")

tmp = tax %>% filter(class == "Lepidosauria")
tmp[, 1:10]
#tax.stat
tax.l.s
genus.tax
genus.tax

# tax.old = read_csv("tables/formatted_species_lower.csv")
# tax.old
# tax.cmp = tax.old %>% left_join(tax, by=c("Species"="str"))
# tax.cmp
# tax.cmp.s = tax.cmp %>% group_by(Major_Class, subphylum, superclass, class, order) %>% summarise(num=n())
# tax.cmp.s
#write_excel_csv(tax.cmp.s, "data/homology/homo.species.tax_cmp.csv")

#tax.cmp %>% filter(Major_Class=="Rodentia", order !="Rodentia") %>% pull(Species)

#tax %>% filter(order=="Testudines")

class2group = c("Actinopteri"="Fish", "Chondrichthyes"="Fish", "Cladistia"="Fish",
                "Aves"="Aves", "Lepidosauria"="Reptilia", "Mammalia"="Mammal",
                "Insecta"="Other", "Saccharomycetes"="Other", "Amphibia"="Other", "Ascidiacea"="Other",
                "Chromadorea"="Other", "Myxini"="Fish", "Hyperoartia"="Fish"
                )
class2group.df = class2group %>% enframe("class", "group1")
tax2g = tax %>% left_join(class2group.df)
unclass1 = tax2g %>% filter(is.na(group1))
unclass1 %>% select(class) %>% table()
unclass1.info = unclass1 %>% select(-superkingdom:-phylum) %>% select(-infraclass:-superorder)

order2group = c("Testudines"="Reptilia", "Crocodylia"="Reptilia", "Coelacanthiformes"="Fish")
tax2g = tax2g %>% left_join(order2group %>% enframe("order", "group2"))
tax2g = tax2g %>% mutate(group.class = ifelse(!is.na(group1), group1, group2))


mammal.order = tax2g %>% filter(group.class=="Mammal") %>% group_by(class, order) %>% summarise(num=n()) %>% arrange(desc(num))
#write_excel_csv(mammal.order, "data/homology/homo.species.mammal2group.csv")
mammal = tax2g %>% filter(group.class=="Mammal")

mammal2group = c("Rodentia"="Rodentia", 
                 #"Proboscidea"="Large","Primates"="Large", "Artiodactyla"="Large", "Carnivora"="Large",
                 "Monotremata"="Prototheria",
                 "Diprotodontia"="Metatheria", "Didelphimorphia"="Metatheria", "Dasyuromorphia"="Metatheria"
                 )
mammal = mammal %>% left_join(mammal2group %>% enframe("order", "groupm"))
mammal = mammal %>% mutate(groupm= ifelse(is.na(groupm), "nrEutheria", groupm))
tax2g =tax2g %>% left_join(mammal %>% select(str, groupm))
write_excel_csv(tax2g, "tables/homology.species2group.detail.csv")

species2g = tax2g %>% select(str, phylum,class, order, group.class, groupm)
species2g
species2g = species2g %>% mutate(group=ifelse(group.class == "Mammal", groupm, group.class))
species2g
species2g = species2g %>% mutate(
  subgroup=ifelse(group.class == "Mammal", order, 
                                ifelse(!is.na(class), class, order)))
species2g
sum(table(species2g$subgroup))
species2g %>% filter(is.na(subgroup))
table(species2g$group)
#write_excel_csv(species2g, "data/homology/homo.species2group.csv")
#write_excel_csv(species2g, "tables/homology.species2group.csv")

# inspect of group details ####
other = tax2g %>% filter(group.class =="Other")
other %>% select(str:class)
tax2g %>% filter(phylum=="Arthropoda")
tax2g %>% select(phylum) %>% distinct()
tax2g %>% filter(subphylum=="Tunicata")
species2g %>% filter(group.class =="Mammal") %>% select(groupm) %>% distinct()

# tax.cmp2 = tax.old %>% left_join(species2g, by=c("Species"="str"))
# tax.cmp2 %>% group_by(group, Major_Class) %>% summarise(num=n()) %>% arrange(desc(num))
# tax.cmp2 %>% filter(group=="nrEutheria", Major_Class != "Large_Mammals")
# tax.cmp2 %>% filter(group=="Fish", Major_Class =="Chordata")

fish = tax2g %>% filter(group.class=="Fish")
fish.info = fish %>% select(str, phylum:order)
fish.info
fish.info.g = fish.info %>% select(-str) %>% distinct()
fish.info.g
fish.info.g$class%>% unique()
fish.info %>% filter(class == "Chondrichthyes")
fish.info %>% filter(superclass == "Sarcopterygii")
