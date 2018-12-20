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

capping <- function(v){
  if (v > 1){return (1)}
  else{return (v)}
  }

plot_fit_heatmap <- function(EspTable,WTresibox,graphname){
  textsize    <- 7
  p <-  ggplot() + 
	  geom_tile(data=EspTable, aes(x=resi,y=aa,fill=mapply(capping,fit_no_antibody))) +
	  scale_fill_gradientn(colours=c("white", "white", "yellow", "red"),
                limits= c(0,1),
		values=rescale(c(0, 0.33, 0.66, 1)),
		guide="colorbar",
                breaks=c(0,0.5,1),
                labels=c('0.0','0.5','1.0'),
		na.value="grey") +
	  theme_classic() +
          theme(axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title=element_text(size=textsize,face="bold"),
                axis.line = element_line(colour = 'black', size = 0),
                legend.title=element_blank(),
	        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
          guides(fill = guide_colorbar(ticks = FALSE, barwidth=0.5, barheight=3,
                                       label.theme = element_text(colour="black",size=textsize,angle=0,face='bold'))) +
          geom_rect(data=WTresibox, size=0.3, fill=NA, colour="black",
                    aes(xmin=x-0.5, xmax=x+0.5, ymin=y-0.5, ymax=y+0.5)) +
          xlab(bquote(bold('HA2 residue')))
  ggsave(graphname,p,width=1.9,height=2.2,dpi=600)
  }

aa_levels  <- rev(c('E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','_'))
pos_levels <- c('42','45','46','47','48','49','52','111')
WTresibox  <- read_tsv('doc/HK68_WTheatmap.tsv')
SEspTable  <- read_tsv('result/Resist_S.tsv') %>%
              mutate(resi=factor(as.character(resi),levels=pos_levels)) %>%
              mutate(aa=factor(as.character(aa),levels=aa_levels))
plot_fit_heatmap(SEspTable,WTresibox,'graph/heatmap_fit_single.png')
