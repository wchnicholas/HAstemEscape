#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(data.table)
require(cowplot)

color_by_resi <- function(resi){
  resi <- str_sub(resi,1,-2)
  if (resi=='(HA2)42'){return ('42')}
  else if (resi=='(HA2)45'){return ('45')}
  else {return ('Other')}
  }

size_by_resi <- function(resi){
  resi <- str_sub(resi,1,-2)
  if (resi=='(HA2)42'){return (3)}
  else if (resi=='(HA2)45'){return (3)}
  else {return (1)}
  }

frac_survive_Perth09_FI6v3 <- read_csv('Bloom_data/Perth09_antibody_FI6v3_median.csv') %>%
                                rename(Perth09_FI6v3=mutfracsurvive)
frac_survive_WSN_FI6v3     <- read_csv('Bloom_data/WSN_antibody_FI6v3_median.csv') %>%
                                rename(WSN_FI6v3=mutfracsurvive)
frac_survive_WSN_CR9114    <- read_csv('Bloom_data/WSN_antibody_CR9114_median.csv') %>%
                                rename(WSN_CR9114=mutfracsurvive)
frac_survive <- frac_survive_Perth09_FI6v3 %>% 
                  inner_join(frac_survive_WSN_FI6v3) %>%
                  inner_join(frac_survive_WSN_CR9114) %>%
                  mutate(variant=mapply(paste,site,mutation,sep='')) %>%
                  select(variant,Perth09_FI6v3,WSN_FI6v3,WSN_CR9114) %>%
                  data.table() %>%
                  melt() %>%
                  mutate(variable=factor(variable,levels=rev(c('Perth09_FI6v3','WSN_FI6v3','WSN_CR9114')))) %>%
                  mutate(color=mapply(color_by_resi,variant)) %>%
                  mutate(size=mapply(size_by_resi,variant)) 

textsize <- 7
sizes <- rev(unique(frac_survive$size))
p  <- ggplot(frac_survive,aes(x=variable,y=value,color=color,size=as.factor(size))) +
        geom_point(position="jitter",alpha=0.7) +
        scale_size_manual(values = sizes/5) +
        scale_color_manual(values=c('blue','red','gray40')) +
        theme(plot.title=element_text(size=textsize,face="bold"),
	      axis.title=element_text(size=textsize,face="bold"),
	      axis.text=element_text(size=textsize,face="bold"),
	      axis.text.x=element_text(face="bold",angle=90,hjust=0,vjust=0.5,size=textsize),
	      legend.title=element_blank(),
	      legend.text=element_text(size=textsize,face="bold"),
	      legend.position='none') +
        ylab("fraction surviving") +
        xlab("")
ggsave('graph/frac_surviving_compare.png',p,height=2.35,width=2.2)
