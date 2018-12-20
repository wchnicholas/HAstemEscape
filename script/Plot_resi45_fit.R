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

fits_table <- read_tsv('result/Fit_compare.tsv')
fits_table_resi45 <- fits_table %>%
                       filter(pos==45) %>%
                       select(aa, HK68_fit_no_antibody, WSN_fit_mock, Perth09_fit_mock) %>%
                       rename(WSN=WSN_fit_mock) %>%
                       rename(Perth09=Perth09_fit_mock) %>%
                       rename(HK68=HK68_fit_no_antibody) %>%
                       data.table() %>%
                       melt()
textsize <- 7
colorscale  <- c(brewer.pal(3,"Accent"))
p <-  ggplot(fits_table_resi45, aes(x=variable, y=value, label=aa, color=variable)) + 
        geom_text(fontface = "bold",position=position_jitter(width=0.4), size=2.5) +
        geom_hline(yintercept = 1, lty=2) +
        scale_color_manual(values=colorscale) +
        theme(plot.title=element_text(size=textsize,face="bold"),
              axis.title=element_text(size=textsize,face="bold"),
              axis.text=element_text(size=textsize,face="bold"),
              legend.title=element_blank(),
              legend.text=element_text(size=textsize,face="bold"),
              legend.position='none') +
        ylab("Relative fitness") +
        xlab("")
ggsave('graph/resi45_fit.png',p,height=1.5,width=4.5)
