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

aa_levels   <- rev(c('E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','_'))
resi_levels <- c('Q42','I45','D46','Q47','I48','N49','L52','T111')
mut_levels <- c()
for (resi in resi_levels){
  for (aa in aa_levels){
    mut_levels <- c(mut_levels,paste(resi,aa,sep=''))
    }
  }

capping1 <- function(v){
  if (v > 1){return (1)}
  else{return (v)}
  }

capping5 <- function(v){
  if (v > 5){return (5)}
  else{return (v)}
  }

plot_fit_heatmap2 <- function(Fit_Table, graphname){
  textsize <- 1
  p <-  ggplot() + 
	  geom_tile(data=Fit_Table, aes(x=mut1,y=mut2,fill=mapply(capping1,score))) +
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
		axis.ticks = element_blank(),
		legend.title=element_blank(),
		panel.border = element_rect(colour = "black", fill=NA, size=0)) +
	  guides(fill = guide_colorbar(ticks = FALSE, barwidth=0.5, barheight=3,
				       label.theme = element_text(colour="black",size=textsize,angle=0,face='bold')))
  ggsave(graphname,p,width=3.5,height=3,dpi=600)
  }

plot_esp_heatmap2 <- function(EspTable,graphname){
  textsize    <- 1
  p <-  ggplot() +
          geom_tile(data=EspTable, aes(x=mut1,y=mut2,fill=mapply(capping5,score))) +
          scale_fill_gradientn(colours=c("green", "white", "white", "white", "purple"),
                limits= c(0,5),
                values=rescale(c(0, 0.8, 1, 2, 5)),
                guide="colorbar",
                breaks=c(0,1,2,3,4,5),
                labels=c('0.0','1.0','2.0','3.0','4.0','5.0'),
                na.value="grey") +
          theme_classic() +
          theme(axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title=element_text(size=textsize,face="bold"),
                axis.line = element_line(colour = 'black', size = 0),
		axis.ticks = element_blank(),
                legend.title=element_blank(),
                panel.border = element_rect(colour = "black", fill=NA, size=0)) +
          guides(fill = guide_colorbar(ticks = FALSE, barwidth=0.5, barheight=3,
                                       label.theme = element_text(colour="black",size=textsize,angle=0,face='bold'))) 
  ggsave(graphname,p,width=3.5,height=3,dpi=600)
  }

Fit_Table <- read_tsv('result/Heatmap_D.tsv') %>%
               mutate(mut1=factor(mut1,level=mut_levels)) %>%
               mutate(mut2=factor(mut2,level=rev(mut_levels)))
fit_no_antibody_Table <- Fit_Table %>% mutate(score=fit_no_antibody)
resist_10ug_CR9114_Table <- Fit_Table %>% mutate(score=resist_10ug_CR9114)
resist_2500ng_FI6v3_Table <- Fit_Table %>% mutate(score=resist_2500ng_FI6v3)
plot_fit_heatmap2(fit_no_antibody_Table, 'graph/heatmap_fit_double.png')
plot_esp_heatmap2(resist_10ug_CR9114_Table, 'graph/heatmap_resist_10ug_CR9114_double.png')
plot_esp_heatmap2(resist_2500ng_FI6v3_Table, 'graph/heatmap_resist_2500ng_FI6v3_double.png')
