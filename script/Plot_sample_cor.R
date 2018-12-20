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
require(cowplot)

coloring <- function(mut){
  if (grepl('_',mut)){return ('Nonsense')}
  else{return ('Missense')}
  }

panel.cor <- function(x, y, digits=2, prefix="", cex.cor){
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  test <- cor.test(x,y) 
  text(0.5, 0.5, txt, cex=1.5) 
  }

SEsp <- read_tsv('result/Resist_S.tsv') %>%
          mutate(color=factor(mapply(coloring,mut),levels=c('Nonsense','Missense'))) %>%
          arrange(desc(color)) %>%
          filter(Input_count >= 20) %>%
          filter(mut!='WT') %>%
          filter(fit_no_antibody >= 0.5) %>%
          select(fit_no_antibody,resist_2ug_CR9114,resist_10ug_CR9114,resist_300ng_FI6v3,resist_2500ng_FI6v3)
          #select(fit_no_antibody,fitA1,fitA5,fitF1,fitF5)
DEsp <- read_tsv('result/Resist_D.tsv') %>%
          mutate(color=factor(mapply(coloring,mut),levels=c('Nonsense','Missense'))) %>%
          arrange(desc(color)) %>%
          filter(Input_count >= 20) %>%
          filter(mut!='WT') %>%
          filter(fit_no_antibody >= 0.5) %>%
          select(fit_no_antibody,resist_2ug_CR9114,resist_10ug_CR9114,resist_300ng_FI6v3,resist_2500ng_FI6v3)
          #select(fit_no_antibody,fitA1,fitA5,fitF1,fitF5)
png('graph/QC_antibodies_cor_S.png',res=300,height=1800,width=1800)
pairs(SEsp, cex=0.5, upper.panel=panel.cor)
dev.off()
png('graph/QC_antibodies_cor_D.png',res=300,height=1800,width=1800)
pairs(DEsp, cex=0.5, upper.panel=panel.cor)
dev.off()
