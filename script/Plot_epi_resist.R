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

plot_epi_resist <- function(data,graphname,xlab,ylab){
  textsize <- 7
  c <- signif(cor(log10(data$resist),log10(data$expect),use="complete.obs",method="pearson"),2)
  print ('############################')
  print (title)
  print (paste('Pearson correlation of log10 fitness:', c, sep=' '))
  p <- ggplot(data,aes(x=log10(expect),y=log10(resist))) +
	 geom_point(size=0.1) +
         geom_hline(yintercept = log10(1), color = 'grey', linetype=2) +
         geom_vline(xintercept = log10(1), color = 'grey', linetype=2) +
	 xlab(xlab) +
	 ylab(ylab) +
         ggtitle(title) +
         theme_cowplot(12) +
	 theme(plot.title=element_blank(),
               axis.title=element_text(size=textsize,face="bold"),
	       axis.text=element_text(size=textsize,face="bold"),
	       legend.title=element_blank(),
	       legend.text=element_text(size=textsize,face="bold"),
	       legend.position='none')
  ggsave(graphname,p,height=2.2,width=2.2)
  }

data <- read_tsv('result/Resist_epi.tsv')
CR9114 <- data %>%
            rename(resist=resist_CR9114) %>%
            rename(expect=expect_resist_CR9114)
FI6v3  <- data %>%
            rename(resist=resist_FI6v3) %>%
            rename(expect=expect_resist_FI6v3)
plot_epi_resist(CR9114,'graph/epi_resist_CR9114.png','Log10 expected relative resistance','Log10 relative resistance')
plot_epi_resist(FI6v3,'graph/epi_resist_FI6v3.png','Log10 expected relative resistance','Log10 relative resistance')
