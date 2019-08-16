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

color_by_resi <- function(resi){
  if (resi==42){return ('42')}
  else if (resi==45){return ('45')}
  else {return ('Other')}
  }

size_by_resi <- function(resi){
  if (resi==42){return (3)}
  else if (resi==45){return (3)}
  else {return (1)}
  }

plot_scatter <- function(fit_table,graphname,xlab,ylab){
  textsize <- 7
  sizes <- unique(fit_table$size)
  names(sizes) <- as.character(unique(fit_table$size))
  cor_summary <- cor.test(log10(fit_table$HK68),log10(fit_table$Other))
  pvalue <- cor_summary$p.value
  cor    <- signif(cor(log10(fit_table$HK68),log10(fit_table$Other)),2)
  p  <- ggplot(fit_table,aes(log10(HK68),log10(Other),color=color,size=as.factor(size))) + 
          geom_point(alpha=0.7) +
          scale_size_manual(values = sizes/5) + 
          scale_color_manual(values=c('blue','red','gray40')) +
          theme(plot.title=element_text(size=textsize,face="bold"),
                axis.title=element_text(size=textsize,face="bold"),
                axis.text=element_text(size=textsize,face="bold"),
                legend.title=element_blank(),
                legend.text=element_text(size=textsize,face="bold"),
                legend.position='none') +
          xlab(xlab) +
          ylab(ylab) +
          geom_hline(aes(yintercept=0), color='green', linetype=2,size=0.2) +
          geom_vline(aes(xintercept=0), color='green', linetype=2,size=0.2)
  print (graphname)
  print (paste('P=',pvalue,sep=''))
  print (paste('R=',cor,sep=''))
  ggsave(graphname, p, width=1.8,height=1.8)
  }

FitTable  <- read_tsv('result/Fit_compare.tsv') %>%
                mutate(color=mapply(color_by_resi,pos)) %>%
                mutate(size=mapply(size_by_resi,pos)) %>%
                filter(HK68_fit_no_antibody!=-1)
fit_table_vsWSN <- FitTable %>%
                     select(HK68_fit_no_antibody, WSN_fit_mock, color, size) %>%
                     rename(HK68=HK68_fit_no_antibody) %>%
                     rename(Other=WSN_fit_mock)
fit_table_vsPerth09 <- FitTable %>%
                         select(HK68_fit_no_antibody, Perth09_fit_mock, color, size) %>%
                         rename(HK68=HK68_fit_no_antibody) %>%
                         rename(Other=Perth09_fit_mock)
fit_table_vsPerth09_FI6v3 <- FitTable %>%
			       select(HK68_fit_2500ng_FI6v3, Perth09_fit_FI6v3_15ug, color, size) %>%
			       rename(HK68=HK68_fit_2500ng_FI6v3) %>%
			       rename(Other=Perth09_fit_FI6v3_15ug)
fit_table_vsWSN_FI6v3 <- FitTable %>%
		           select(HK68_fit_2500ng_FI6v3, WSN_fit_FI6v3_200ng, color, size) %>%
		           rename(HK68=HK68_fit_2500ng_FI6v3) %>%
		           rename(Other=WSN_fit_FI6v3_200ng)
fit_table_vsWSN_CR9114 <- FitTable %>%
		            select(HK68_fit_10ug_CR9114, WSN_fit_CR9114_100ng, color, size) %>%
		            rename(HK68=HK68_fit_10ug_CR9114) %>%
		            rename(Other=WSN_fit_CR9114_100ng)
plot_scatter(fit_table_vsWSN,'graph/scatter_fit_WSN_compare.png','HK68 (no antibody)','WSN (no antibody)')
plot_scatter(fit_table_vsPerth09,'graph/scatter_fit_Perth09_compare.png','HK68 (no antibody)','Perth09 (no antibody)')
plot_scatter(fit_table_vsPerth09_FI6v3,'graph/scatter_fit_Perth09_FI6v3_compare.png','HK68 (2.5 ug/mL FI6v3)','Perth09 (15 ug/mL FI6v3)')
plot_scatter(fit_table_vsWSN_FI6v3,'graph/scatter_fit_WSN_FI6v3_compare.png','HK68 (2.5 ug/mL FI6v3)','WSN (0.2 ug/mL FI6v3)')
plot_scatter(fit_table_vsWSN_CR9114,'graph/scatter_fit_WSN_CR9114_compare.png','HK68 (10 ug/mL CR9114)','WSN (0.1 ug/mL CR9114)')
