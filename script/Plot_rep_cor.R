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

Plot_RepvsRep <- function(lib,xlab,ylab,title){
  textsize <- 7
  c <- signif(cor(log10(lib$Rep1fit),log10(lib$Rep2fit),use="complete.obs",method="pearson"),2)
  print ('############################')
  print (title)
  print (paste('Pearson correlation of log10 fitness:', c, sep=' '))
  p <- ggplot(lib,aes(x=log10(Rep1fit),y=log10(Rep2fit),color=color)) +
	 geom_point(size=0.1) +
         geom_hline(yintercept = log10(1), color = 'grey', linetype=2) +
         geom_vline(xintercept = log10(1), color = 'grey', linetype=2) +
	 xlab(xlab) +
	 ylab(ylab) +
	 scale_color_manual(values=c('#A0A000AF','#7001BAAF')) +
         ggtitle(title) +
	 theme(plot.title=element_blank(),
               axis.title=element_text(size=textsize,face="bold"),
	       axis.text=element_text(size=textsize,face="bold"),
	       legend.title=element_blank(),
	       legend.text=element_text(size=textsize,face="bold"),
	       legend.position='none') +
         scale_x_continuous(breaks = c(-2,-1,0,1,2))
  return (p)
  }

Slib <- read_tsv('result/Fitness_S.tsv') %>%
          mutate(color=factor(mapply(coloring,mut),levels=c('Nonsense','Missense'))) %>%
          arrange(desc(color)) %>%
          filter(Input_count >= 20)
Dlib <- read_tsv('result/Fitness_D.tsv') %>%
          mutate(color=factor(mapply(coloring,mut),levels=c('Nonsense','Missense'))) %>%
          arrange(desc(color)) %>%
          filter(Input_count >= 20)
Slib_A0 <- Slib %>% 
             rename(Rep1fit=rep1_no_antibody_fit) %>%
             rename(Rep2fit=rep2_no_antibody_fit)
Slib_A1 <- Slib %>% 
             rename(Rep1fit=rep1_2ug_CR9114_fit) %>%
             rename(Rep2fit=rep2_2ug_CR9114_fit)
Slib_A5 <- Slib %>% 
             rename(Rep1fit=rep1_10ug_CR9114_fit) %>%
             rename(Rep2fit=rep2_10ug_CR9114_fit)
Slib_F1 <- Slib %>% 
             rename(Rep1fit=rep1_300ng_FI6v3_fit) %>%
             rename(Rep2fit=rep2_300ng_FI6v3_fit)
Slib_F5 <- Slib %>% 
             rename(Rep1fit=rep1_2500ng_FI6v3_fit) %>%
             rename(Rep2fit=rep2_2500ng_FI6v3_fit)
Dlib_A0 <- Dlib %>% 
             rename(Rep1fit=rep1_no_antibody_fit) %>%
             rename(Rep2fit=rep2_no_antibody_fit)
Dlib_A1 <- Dlib %>% 
             rename(Rep1fit=rep1_2ug_CR9114_fit) %>%
             rename(Rep2fit=rep2_2ug_CR9114_fit)
Dlib_A5 <- Dlib %>% 
             rename(Rep1fit=rep1_10ug_CR9114_fit) %>%
             rename(Rep2fit=rep2_10ug_CR9114_fit)
Dlib_F1 <- Dlib %>% 
             rename(Rep1fit=rep1_300ng_FI6v3_fit) %>%
             rename(Rep2fit=rep2_300ng_FI6v3_fit)
Dlib_F5 <- Dlib %>% 
             rename(Rep1fit=rep1_2500ng_FI6v3_fit) %>%
             rename(Rep2fit=rep2_2500ng_FI6v3_fit)
p1  <- Plot_RepvsRep(Slib_A0,'Replicate 1','Replicate 2',expression(bold('no antibody')))
p2  <- Plot_RepvsRep(Slib_A1,'Replicate 1','Replicate 2',expression(bold("CR9114 (2"~"\u03bc"*"g/mL)")))
p3  <- Plot_RepvsRep(Slib_A5,'Replicate 1','Replicate 2',expression(bold("CR9114 (10"~"\u03bc"*"g/mL)")))
p4  <- Plot_RepvsRep(Slib_F1,'Replicate 1','Replicate 2',expression(bold("FI6v3 (0.3"~"\u03bc"*"g/mL)")))
p5  <- Plot_RepvsRep(Slib_F5,'Replicate 1','Replicate 2',expression(bold("FI6v3 (2.5"~"\u03bc"*"g/mL)")))
p6  <- Plot_RepvsRep(Dlib_A0,'Replicate 1','Replicate 2',expression(bold('no antibody')))
p7  <- Plot_RepvsRep(Dlib_A1,'Replicate 1','Replicate 2',expression(bold("CR9114 (2"~"\u03bc"*"g/mL)")))
p8  <- Plot_RepvsRep(Dlib_A5,'Replicate 1','Replicate 2',expression(bold("CR9114 (10"~"\u03bc"*"g/mL)")))
p9  <- Plot_RepvsRep(Dlib_F1,'Replicate 1','Replicate 2',expression(bold("FI6v3 (0.3"~"\u03bc"*"g/mL)")))
p10 <- Plot_RepvsRep(Dlib_F5,'Replicate 1','Replicate 2',expression(bold("FI6v3 (2.5"~"\u03bc"*"g/mL)")))
p <- grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,nrow=2,ncol=5)
ggsave('graph/QC_replicates_cor.png',p,height=2.8,width=6.5)
