##This is based (actually almost verbatim atm) on this dudes workflow: https://www.lesleysitter.com/2019/08/23/roary-analysis/ 
## I have fucked with that dudes pipeline like a lot at this point, but it still is ocasionally helpful to go read it if you're stuck\

#### The whole point of this is to translate the presence absence csv from roary into a visual output, ordered to PAO1 to compare genomes and hopefully find deletions. 


table_input <- read.csv("rorarydfast .csv", sep=",", na.strings = c("","NA"))
table_input_df <- as.data.frame(table_input)

##when the csv is ordered, it does not retain colnames, this is me fixing that
##fixed in excel stage legacy code 
##table_input_ordered <- read.csv("4iso - 4iso.csv", sep=",", na.strings = c("","NA"))
##table_input_ordered_df <- as.data.frame(table_input_ordered)
##colnames(table_input_ordered_df) <- colnames(table_input_df)
##View(table_input_ordered_df)

## removes all the extra stuff that we aren't interested in for this visualisation
tableValues <- within(table_input, rm("Annotation", "Non.unique.Gene.name","No..isolates","No..sequences","Avg.sequences.per.isolate","Genome.Fragment","Order.within.Fragment","Accessory.Fragment","Accessory.Order.with.Fragment","QC","Min.group.size.nuc","Max.group.size.nuc","Avg.group.size.nuc"))
##View(tableValues)

## This turns the unhelpful gene ID's that roary prints into 1 for presence and 0 for absence 
abscence_presence <- as.matrix(tableValues[,-1])
rownames(abscence_presence) <- tableValues[,1]
abscence_presence[is.na(abscence_presence)] <- 0
abscence_presence[which(abscence_presence!=0)] <- 1

##makes the matrix into a useable object with useful roawnames and colnames
a_p_matrix <- mapply(abscence_presence, FUN=as.numeric)
a_p_matrix <- matrix(data=a_p_matrix, ncol=length(colnames(abscence_presence)), nrow=length(row.names(abscence_presence)))
row.names(a_p_matrix) <- row.names(abscence_presence)
colnames(a_p_matrix) <- colnames(abscence_presence)

##This is to tidy up long gross file names when we are at the heatmap stage. Make sure you have the same length :) 
##colnames(a_p_matrix) <- c("020MIC", "13121.1", "15108", "1709.12", "2192", "39016", "39177", "57P31PA", "5C1_S24", "679", "968333S", "A5803", "AA2", "AES.1R_2", "AMT0023.30", "AMT0023.34", "AMT0060.1", "AMT0060.2", "AMT0060.3", "AUS058", "AUS066", "AUS088", "AUS089", "AUS23", "AUS52", "C3179", "CHA", "CPHL9433", "DK2", "DUN.001C", "DUN.003B", "DUN.004", "ILPAO1", "IST27", "IST27N", "Jpn1563", "KK1", "LES400", "LES431", "LESB58", "LMG14084", "Mi162_2", "NH57388A","PAO1", "PA7", "PAK", "Pr335", "TBCF10839", "U018a", "PA14", "PAO1_roary" )
##View(a_p_matrix)
##View(abscence_presence)

##to remove genes that aren't present in PAO1, I will turn my nice useable matrix back into a dataframe, subset it and then turn it back into a matrix and it should retain its useability 
cheating <-as.data.frame(a_p_matrix)
##View(cheating)
##Removes any genes that only are present outside of PAO1
subsetcheating1 <- cheating[!(cheating$DFAST_PA01.fasta==0),]
##View(subsetcheating1)
##make this back into a matrix
altered_ap_matrix <- as.matrix(subsetcheating1)
##View(altered_ap_matrix)

##I was removing genes based on prevelance, but in terms of finding deletions that have relatively similar genome coordinates as will be presented by this visualisation I am stopping doing that and instead just using the above. Code kept below in case I want to go back to doing it the way I was. 
##Removes genes that are in every strain
##subsercheating2 <- subsetcheating1[!(subsetcheating1$newCol >= 143),]
#View(subsercheating2)
##Removes my cheating column
##altered_ap <- subsercheating2[,1:143]
#altered_ap <- as.numeric(altered_ap)
#View(altered_ap)

##heatmap time. Importantly when you export make it like 1000 in by 500 in to try and avoid axis smooshing

heatmap(altered_ap_matrix, keep.dendro = TRUE, hclustfun = hclust, verbose = TRUE, Rowv = NA, cexCol = 0.5, col = c("#FFE986","#FF736E"), margins = c(4,0.5))

##This second half I have not touched at all pretty much. Certainly could though! but ignore it, its here till I decide what I want to do with it.
genomes_count <- length(colnames(a_p_matrix))

abscence_presence <- cbind(a_p_matrix, rowSums(a_p_matrix))

summary_table <- matrix(data=NA, nrow=3, ncol=length(colnames(abscence_presence)))
colnames(summary_table) <- colnames(abscence_presence)
rownames(summary_table) <- c("Total_genes","Unique_genes","Core_genes")

summary_table[1,] <- colSums(abscence_presence)
summary_table[2,] <- colSums(abscence_presence[which(abscence_presence[,ncol(abscence_presence)] == 1),])
summary_table[3,] <- colSums(abscence_presence[which(abscence_presence[,ncol(abscence_presence)] >= (genomes_count*0.95)),])
summary_table <- summary_table[,-ncol(summary_table)]

average_table <- data.frame(x=1:6, y=1:6, z=1:6)

average_table[,1] <- c("Total genes analyzed","Orthologous groups","Average gene count","Average core genes","Average unique genes","Total unique genes")
average_table[1,2] <- sum(summary_table[1,])
average_table[2,2] <- length(rownames(abscence_presence))
average_table[3,2] <- median(summary_table[1,])
average_table[4,2] <- length(rownames(abscence_presence[which(abscence_presence[,ncol(abscence_presence)] >= (genomes_count*0.95)),]))
average_table[5,2] <- round(length(rownames(abscence_presence[which(abscence_presence[,ncol(abscence_presence)] == 1),]))/length(colnames(abscence_presence)))
average_table[6,2] <- length(rownames(abscence_presence[which(abscence_presence[,ncol(abscence_presence)] == 1),]))


library("reshape2")
melt_summary_table <- melt(summary_table)
melt_summary_table <- melt_summary_table[order(melt_summary_table$value),]

library(tidyverse)

p1 <- ggplot(melt_summary_table, aes(x = reorder(Var2, value), y = value)) + 
  geom_bar( stat = 'identity') + 
  facet_grid(. ~ Var1, scales = "free_x") + 
  xlab("Genomes") +
  ylab("Count") +
  coord_flip()
p1

p2 <- ggplot(data=average_table[-c(1,2,6),], aes(x=x, y=y))+
  geom_bar(stat = 'identity') +
  theme (axis.text.x=element_text(angle=90,hjust=1,vjust=0.3),
         axis.title.x = element_blank()) +
  geom_text(aes(y = 10, label = paste("N =" ,y),vjust = 0), colour = "white", show.legend=FALSE) +
  ylab("Count")
p2

library(grid)
library(gridExtra)
t1  <- textGrob(paste(c("Total number of genomes:\n",
                        length(colnames(summary_table)),
                        "\n\nNumber of analyzed genes:\n",
                        as.numeric(average_table[1,2]),
                        "\n\nTotal orthologous groups\n",
                        as.numeric(average_table[2,2]),
                        "\n\nTotal unique genes\n",
                        as.numeric(average_table[6,2])), collapse = " ")) 
t1

lay <- rbind(c(1,1,2),
             c(1,1,2),
             c(1,1,3),
             c(1,1,3))

grid.arrange(p1,p2,t1, layout_matrix = lay)
