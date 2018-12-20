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
  if (v > 5){return (5)}
  else{return (v)}
  }

extract_SEscape <- function(DEspTable, mut_bg, aa_levels, pos_levels, SEspTable){
  pos_bg <- str_sub(mut_bg,2,-2)
  new_pos_levels <- c()
  for (pos in pos_levels){
    if (pos!=pos_bg){
      new_pos_levels <- c(new_pos_levels, pos)
      }
    }
  A0 <- filter(SEspTable, mut==mut_bg)$fit_no_antibody
  A5 <- filter(SEspTable, mut==mut_bg)$resist_10ug_CR9114
  F5 <- filter(SEspTable, mut==mut_bg)$resist_2500ng_FI6v3
  A0_list   <- rep(A0,8)
  A5_list   <- rep(A5,8)
  F5_list   <- rep(F5,8)
  mut_list  <- c('Q42Q','I45I','D46D','Q47Q','I48I','N49N','L52L','T111T')
  WT_table  <- data.frame(cbind(mut_list,A0_list,A5_list,F5_list)) %>%
                 rename(mut=mut_list) %>%
                 rename(fit_no_antibody=A0_list) %>%
                 rename(resist_10ug_CR9114=A5_list) %>%
                 rename(resist_2500ng_FI6v3=F5_list) %>%
                 filter(str_sub(mut,2,-2)!=str_sub(mut_bg,2,-2))
  resi_SEspTable_1 <- DEspTable %>% 
                      filter(mut1==mut_bg) %>% 
                      rename(mut=mut2) %>%
                      select(-mut1)
  resi_SEspTable_2 <- DEspTable %>%
                      filter(mut2==mut_bg) %>%
                      rename(mut=mut1) %>%
                      select(-mut2)
  resi_SEspTable   <- rbind(resi_SEspTable_1, resi_SEspTable_2, WT_table) %>%
                      mutate(resi=str_sub(mut,2,-2)) %>% 
                      mutate(aa=str_sub(mut,-1)) %>%
                      mutate(resi=factor(as.character(resi),levels=new_pos_levels)) %>%
                      mutate(aa=factor(as.character(aa),levels=aa_levels)) %>%
                      mutate(fit_no_antibody=mapply(as.numeric,fit_no_antibody)) %>%
                      mutate(resist_10ug_CR9114=mapply(as.numeric,resist_10ug_CR9114)) %>%
                      mutate(resist_2500ng_FI6v3=mapply(as.numeric,resist_2500ng_FI6v3))
  return (resi_SEspTable)
  }

plot_fit_heatmap <- function(EspTable,pos_bg, WTresibox,graphname){
  textsize    <- 7
  WTresibox   <- WTresibox %>%
                   filter(str_sub(resi,2,-1)!=pos_bg) %>%
                   mutate(x=seq(1,length(.$x)))
  if (pos_bg != ''){
    WTresibox <- WTresibox %>%
                   mutate(y=y-1)
    }
  p <-  ggplot() + 
	  geom_tile(data=EspTable, aes(x=resi,y=aa,fill=mapply(capping,X))) +
	  scale_fill_gradientn(colours=c("green", "white", "white", "white", "purple"),
                limits= c(0,5),
		values=rescale(c(0, 0.5, 1, 2, 5)),
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
EspF5_table <- SEspTable %>%
                 mutate(X=resist_2500ng_FI6v3)
EspA5_table <- SEspTable %>%
                 mutate(X=resist_10ug_CR9114)
DEspTable  <- read_tsv('result/Heatmap_D.tsv')
EspTable_I45F <- extract_SEscape(DEspTable, 'I45F', aa_levels, pos_levels, SEspTable)
EspA5_I45F    <- EspTable_I45F %>% mutate(X=resist_10ug_CR9114)                   
EspF5_I45F    <- EspTable_I45F %>% mutate(X=resist_2500ng_FI6v3)                   
EspTable_I45Y <- extract_SEscape(DEspTable, 'I45Y', aa_levels, pos_levels, SEspTable)
EspA5_I45Y    <- EspTable_I45Y %>% mutate(X=resist_10ug_CR9114)                   
EspF5_I45Y    <- EspTable_I45Y %>% mutate(X=resist_2500ng_FI6v3)                   
EspTable_I45W <- extract_SEscape(DEspTable, 'I45W', aa_levels, pos_levels, SEspTable)
EspA5_I45W    <- EspTable_I45W %>% mutate(X=resist_10ug_CR9114)                   
EspF5_I45W    <- EspTable_I45W %>% mutate(X=resist_2500ng_FI6v3)                   
plot_fit_heatmap(EspF5_table,'',WTresibox,'graph/heatmap_FI6v3esp_single.png')
plot_fit_heatmap(EspA5_table,'',WTresibox,'graph/heatmap_CR9114esp_single.png')
plot_fit_heatmap(EspF5_I45F,'45',WTresibox,'graph/heatmap_FI6v3esp_single_onI45F.png')
plot_fit_heatmap(EspA5_I45F,'45',WTresibox,'graph/heatmap_CR9114esp_single_onI45F.png')
plot_fit_heatmap(EspF5_I45Y,'45',WTresibox,'graph/heatmap_FI6v3esp_single_onI45Y.png')
plot_fit_heatmap(EspA5_I45Y,'45',WTresibox,'graph/heatmap_CR9114esp_single_onI45Y.png')
plot_fit_heatmap(EspA5_I45W,'45',WTresibox,'graph/heatmap_CR9114esp_single_onI45W.png')
plot_fit_heatmap(EspF5_I45W,'45',WTresibox,'graph/heatmap_FI6v3esp_single_onI45W.png')
