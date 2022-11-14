#####################################################################
####
####  Timing of Gains with Chemotherapy Barcoding - 
####
#### Creator: Benjamin Diamond and Francesco Maura
####
#####################################################################

#### R code to identify mutational signatures in mutations that are duplicated (or nonduplicated) 
#### within large chromosomal gains (>1 Mb)

library(mmsig)
library(dplyr)
library(ggplot2)
library(stringr)

#####load in sample copy number data and phylogeny data (e.g. Pyclone output)
##### this is required to select the clonal variants and seprated them from subclonal

setwd("~/dupli.sampledata")

all_dp2<- read.table("dupli.phylo.sample.txt") # clonality data generated with pyclone
batt3<- read.table("dupli.cnv.sample.txt") # copy number generated with battenberg

# intersect mutations and copy number and select duplicated copy number segments, 
# assign duplicated annotation per corrected VAF (VAF/purity)

all_dupli<- batt3 %>% dplyr::left_join(all_dp2, by = c("sample", "chr")) %>%
  dplyr::filter(pos<endpos, pos>startpos) %>% 
  dplyr::filter(tot>2.5 | (tot>1.8 & min==0)) %>%
  dplyr::mutate(dupli = ifelse(vaf_crr>0.5, "duplicated",
                               ifelse(vaf_crr<=0.2, "subclonal", "nonduplicated")))

#remove the short segment in this case

all_dupli <- all_dupli %>% dplyr::filter(chr != "4")


# to collapse multiple gains for the purposes of signature analysis, they should be estimated to have occurred within the same time window
# we use molecular time for this purpose as published previously in multiple studies including 
# Rustad et al. Nat Comm. 2020;Maura et al. Nat Comm, 2019; Oben et al. Nat Comm 2021
# the full code is avilable on our lab github: # https://github.com/UM-Myeloma-Genomics/mol_time
# output is included in sample data for reference to this case 

gg_th <- (panel_theme = theme_bw() + theme(
  text = element_text(color = "black"),
  axis.text = element_text(color = "black"),
  panel.border = element_blank(),
  legend.position = 'none',
  legend.key.size = unit(3, 'mm'),
  legend.title = element_blank(),
  legend.direction = 'horizontal',
  legend.text = element_text(size = 8),
  plot.subtitle = element_text(hjust = 0.5, size = 8),
  plot.title = element_text(size = 10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  axis.text.y = element_text(size = 7),
  axis.text.x = element_text(size = 7),
  axis.title = element_text(size = 8),
  axis.line = element_line(),
  plot.margin = unit(c(3,2,2,2), 'pt'))) 

### this defines which chromsomal gains are acquired in the same time window, and which is the earliest time window. 

tmnmt<-read.delim("IID_H130563_T05_01_WG01_cluster_summary.txt", stringsAsFactors = F) 
tmnmt$label<-paste(tmnmt$code, tmnmt$chr,sep="_")
tmnmta<- tmnmt %>%
  dplyr::select(3,4,5,15)
tmnmtb<-tmnmt %>%
  dplyr::select(6,7,8,15) %>%
  dplyr::rename(plot_dot1=plot_dot2, IC_down_cn1=IC_down_cn2, IC_up_cn1=IC_up_cn2)


tmnmt2<- rbind(tmnmta,tmnmtb)[1:5,]
tmnmt2$col<-c("a","a","a","a","b")

ggplot(data=tmnmt2[c(1:4),], aes(x=label, y=plot_dot1, fill=col)) +  geom_dotplot(binaxis="y", stackdir="center", dotsize = 1, position = position_dodge(0.5), stackratio = 1) +
  scale_fill_brewer(palette = "Dark2") + geom_errorbar(aes(ymin=IC_down_cn1, ymax=IC_up_cn1), size = 0.5, width=0.3) + gg_th +
  labs(title = element_blank()) + scale_y_continuous(limits = c(0,1, breaks = seq(0, 1, by = 0.2))) + xlab(label = element_blank()) + ylab(label = "molecular time")

### mutations acquired in gains occurred in same time window are collpased together


############################################################
#####
##### Duplicated/Non-duplicated chemo mutations
#####
###########################################################

#prep data for mmSig
all_dupli$sample<- paste(all_dupli$sample, all_dupli$dupli, sep = " ")
all_dupli$ref<- str_split_fixed(all_dupli$mutation_id, "_", 9)[,8]
all_dupli$alt<- str_split_fixed(all_dupli$mutation_id, "_", 9)[,9]
all_dupli<- all_dupli[,c(1,2,18,24,25)]
colnames(all_dupli)<- c("sample","chr", "pos", "ref", "alt")
all_dupli$chr<- paste("chr", all_dupli$chr, sep = "")


# create per-sample signature catalogue for signatures known to be active in the tumor - 
# here we use our prior global mmsig output (see extended data table 4 for rest of samples)
#
# mmsig original code can be found: https://github.com/UM-Myeloma-Genomics/mmsig

comb<- read.delim("dupli.sig.catalogue.txt", sep = "\t", stringsAsFactors = F)
colnames(comb)<- gsub("[.]", "-", colnames(comb))  

#generate the catalogue
sig.list.dupli<- list()
for(i in (1:nrow(comb))){
  
  signature_pos<- colnames(comb)[7:16][which(comb[i,7:16]!=0)]
  signature_pos[signature_pos ==  "AML-HSC"]<-"SBS5" ##if running on myeloid samples, rename HSC signature as SBS5 for mmsig to accept as clock-like
  sig.list.dupli[[paste(comb$sample[i], "duplicated", sep=" ")]]<- c(paste(comb$sample[i], "duplicated", sep=" ") , signature_pos)
  sig.list.dupli[[paste(comb$sample[i], "nonduplicated", sep=" ")]]<- c(paste(comb$sample[i], "nonduplicated", sep=" ") , signature_pos)
  sig.list.dupli[[paste(comb$sample[i], "subclonal", sep=" ")]]<- c(paste(comb$sample[i], "subclonal", sep=" ") , signature_pos)
  
}

#load reference signatures
sig_ref1<- read.table("sig_ref.txt",stringsAsFactors = F, sep = "\t")
colnames(sig_ref1)<- gsub("[.]", "-", colnames(sig_ref1))   

#run mmsig for mutational signatures estimation
sig_out_dupli <- mm_fit_signatures(muts.input=all_dupli, 
                                   sig.input=sig_ref1, 
                                   input.format = "vcf",
                                   sample.sigt.profs = sig.list.dupli, 
                                   strandbias = TRUE,
                                   bootstrap = TRUE,
                                   iterations = 10, # 1000 iterations recommended for stable results
                                   dbg=FALSE) 



colnames(sig_out_dupli$estimate)[6]<- "SBS-HSC" #revert back to original signature name
sig_out_dupli$estimate

# In this example, we see a melphalan signature (SBS-MM1) among duplicated mutations within two CN-LOH. 
# This means that chemotherapy exposure must have occurred prior to these gains for this to be possible


# We can plot the cancer cell fraction of each mutation within each gain for visualization of duplicated and non duplicated mutations 
# From the output of molecular time analysis

f<- list.files()[grep("cluster_muts.txt", list.files())]
all_mutccf<- list()
for(i in (1:length(f))){
  
  file_int<- read.delim(f[i], stringsAsFactors = F)
  file_int$sample<- gsub("_cluster_muts.txt","", f[i])
  all_mutccf[[i]]<- file_int
}
all_mutccf<- do.call(rbind, all_mutccf)
all_mutccf$code<- as.character(all_mutccf$code)
all_mutccf<- all_mutccf %>% dplyr::mutate(code = ifelse(code == "2", "Duplicated", 
                                                        ifelse(code == "1", "Non-duplicated",
                                                               ifelse(code == "3", "NA", code))))
annot_point = list(c("Duplicated" = "#FBB4AE", "NA" = "#DDDDDD", "Non-duplicated" = "#BEAED4")) 

gg_th <- (panel_theme = theme_bw() + theme(
  panel.border = element_blank(),
  text = element_text(color = "black"),
  axis.text = element_text(color = "black"),
  legend.position = 'right',
  legend.key.size = unit(3, 'mm'),
  legend.direction = 'horizontal',
  legend.title = element_blank(),
  legend.text = element_text(size = 7), 
  plot.subtitle = element_text(hjust = 0.5, size = 8),
  plot.title = element_text(face = 'bold', hjust = 0, vjust = -2, size = 12),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  axis.text.y = element_text(size = 8),
  axis.text.x = element_text(size = 7), 
  axis.title = element_text(size = 8),
  axis.line = element_line(),
  plot.margin = unit(c(3,2,2,2), 'pt'))) 



ggplot(all_mutccf[all_mutccf$sample == "IID_H130563_T05_01_WG01_5",], aes(x = Pos, y = ccf)) + geom_point(aes(color = code)) + gg_th +
  scale_color_manual(values =c("Duplicated" = "#FBB4AE", "NA" = "#DDDDDD", "Non-duplicated" = "#BEAED4")) + 
  scale_y_continuous(limits=c(0, 100), breaks = seq(0,100, by = 10)) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  xlab(label = "position")


# we can plot the 96-class profile of mutations within the gains as well 
profiles_prepost<- sig_out_dupli$mutmatrix
par(mfrow=c(1,1), mar=c(8,5,10,1))
mm_col <- c(rep("lightskyblue",16), rep("black",16), rep("firebrick2",16),rep("gray88",16),rep("darkolivegreen3",16),
            rep("lightpink1",16))
jj<- as.numeric(profiles_prepost[,1])
barplot(jj, col=mm_col, names.arg = rownames(profiles_prepost), las=2, border=F, cex.axis=1.5)