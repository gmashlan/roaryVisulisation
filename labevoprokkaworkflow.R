##This script assumes you have an input file from roary. (absence presence csv) after reordering contigs with mauve, annotating using DFAST/prokka and then running roary. The output of this script is an ordered relative to reference geneome heatmap-esq plot, which can be zoomed in on/manipulated to serve your purposes. Props to this workflow:  https://www.lesleysitter.com/2019/08/23/roary-analysis/ which this script is based on.
library(tidyverse)
table_input_labevoP <- read.csv("prokkalabevoPAO1.csv", sep=",", na.strings = c("","NA"))
table_input_df_labevoP <- as.data.frame(table_input_labevoP)
## removes all the extra stuff that we aren't interested in for this visualisation
tableValues_labevoP <- within(table_input_labevoP, rm("Annotation", "Non.unique.Gene.name","No..isolates","No..sequences","Avg.sequences.per.isolate","Genome.Fragment","Order.within.Fragment","Accessory.Fragment","Accessory.Order.with.Fragment","QC","Min.group.size.nuc","Max.group.size.nuc","Avg.group.size.nuc"))
##Orders your genes relative to reference genome
orderedvalues_labevoP <- tableValues_labevoP[order(tableValues_labevoP$Pseudomonas_aeruginosa_PAO1_107.fna),] ##need to change your col name to whatever your ref col is
## This turns the unhelpful gene ID's that roary prints into 1 for presence and 0 for absence 
abscence_presence_labevoP <- as.matrix(orderedvalues_labevoP[,-1])
rownames(abscence_presence_labevoP) <- orderedvalues_labevoP[,1]
abscence_presence_labevoP[is.na(abscence_presence_labevoP)] <- 0
abscence_presence_labevoP[which(abscence_presence_labevoP!=0)] <- 1
##makes the matrix into a useable object with useful rownames and colnames
a_p_matrix_labevoP <- mapply(abscence_presence_labevoP, FUN=as.numeric)
a_p_matrix_labevoP <- matrix(data=a_p_matrix_labevoP, ncol=length(colnames(abscence_presence_labevoP)), nrow=length(row.names(abscence_presence_labevoP)))
row.names(a_p_matrix_labevoP) <- row.names(abscence_presence_labevoP)
colnames(a_p_matrix_labevoP) <- colnames(abscence_presence_labevoP)
##This is to tidy up long gross file names when we are at the heatmap stage. Make sure you have the same length :)
colnames(a_p_matrix_labevoP) <- c("Aus088", "AUS089", "AUS52", "DUN003B", "DUN004", "ILPAO1", "MC_post_del", "MC_pre_del", "MF_post_del", "MF_pre_del", "MH_post_del", "MH_pre_del", "MJ_post_del", "MJ_pre_del", "MK_post_del", "MK_pre_del", "M_post_del", "M_pre_del", "PAO1")
#to remove genes that aren't present in PAO1, I will turn my nice useable matrix back into a dataframe, subset it and then turn it back into a matrix and it should retain its useability
cheating_labevoP <-as.data.frame(a_p_matrix_labevoP)
##Removes any genes that only are present outside of PAO1/ref
subsetcheating1_labevoP <- cheating_labevoP[!(cheating_labevoP$PAO1==0),]##need to change this to whatever ref is too
subsetcheating1_labevoP <- subsetcheating1_labevoP %>% relocate(PAO1, .before =Aus088) ## this is to re-order so your ref is first ##Turn it back into a matrix
altered_ap_matrix_labevoP <- as.matrix(subsetcheating1_labevoP)
##make a really basic gene rank % lost basically
generanks_labevoP <-as.data.frame(rowSums(subsetcheating1_labevoP))/(ncol(subsetcheating1_labevoP))
generanks_labevoP$gene <- rownames(subsetcheating1_labevoP)
colnames(generanks_labevoP)<- c("rank", "gene")
##remove all 1s/present in most (cut off is 90% atm)
interestinggeneranks_labevoP <- generanks_labevoP[!(generanks_labevoP$rank)>=0.90,]
interestinggeneranks_labevoP <- as.data.frame(interestinggeneranks_labevoP)
orderedgeneranks_labevoP <- interestinggeneranks_labevoP[order(interestinggeneranks_labevoP$rank, decreasing = FALSE),]
first100ranks_labevoP <- orderedgeneranks_labevoP[1:100,]
##plot your generanks
pdf(file = "first100rankslabevoProkka.pdf", height = 10)
ggplot(first100ranks_labevoP, aes(x=gene, y= rank)) + geom_bar(stat = "identity", width = 0.6) +coord_flip()
dev.off()
##Plot your absence/presence big online resource version
##data needs to be long for ggplot to handle
melted_matrix_labevoP <- reshape2::melt(altered_ap_matrix_labevoP)
pdf(file = "prokkalabevo_big.pdf", height=100, width = 10)
ggplot(melted_matrix_labevoP, aes(x= Var2, y=Var1, fill=c("dark green", "grey")[value+1]))+
  geom_raster()+
  scale_fill_identity()+
  theme(axis.text.x=element_text(size = 4, angle=90),
        axis.text.y = element_text(size = 1.5), ##this needs to be made into element_blank if not printing a massive figure
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) 
dev.off()
##plot your absence presence compressed overview version
pdf(file = "prokkalabevo_small.pdf", height=10, width = 10)
ggplot(melted_matrix_labevoP, aes(x= Var2, y=Var1, fill=c("dark green", "grey")[value+1]))+
  geom_raster()+
  scale_fill_identity()+
  theme(axis.text.x=element_text(size = 4, angle=90),
        axis.text.y = element_blank(), ##this needs to be made into element_blank if not printing a massive figure
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) 
dev.off()
###plot your relavent zoomed location 
##these are how you find the numbers to subset based on, look at the *_big output to find the genes you want. 
which(orderedvalues_labevoP$Gene == 'group_3230')  
which (orderedvalues_labevoP$Gene == 'group_5882')
zoomed_ap_cheating_labevoP <- subsetcheating1_labevoP[1829:2464,] ##put the outputs of the which in here :) 
##back into a matrix
zoomed_Ap_matrix_labevoP <- as.matrix(zoomed_ap_cheating_labevoP)
zoomed_melted_matrix_labevoP <- reshape2::melt(zoomed_Ap_matrix_labevoP)
##plot time
pdf(file = "zoomedprokkalabevo_small.pdf")
ggplot(zoomed_melted_matrix_labevoP, aes(x= Var2, y=Var1, fill=c("dark green", "grey")[value+1]))+
  geom_raster()+
  scale_fill_identity()+
  theme(axis.text.x=element_text(size = 10, angle=90),
        axis.text.y = element_blank(), ##this needs to be made into element_blank if not printing a massive figure
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())
dev.off()