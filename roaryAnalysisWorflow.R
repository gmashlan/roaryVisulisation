##This is based (actually almost verbatim atm) on this dudes workflow: https://www.lesleysitter.com/2019/08/23/roary-analysis/ 
## I have fucked with that dudes pipeline like a lot at this point, but it still is ocasionally helpful to go read it if you're stuck


table_input <- read.csv("DFASTandMauvePresAbs.csv", sep=",", na.strings = c("","NA"))
table_input_df <- as.data.frame(table_input)
##when the csv is ordered, it does not retain colnames, this is me fixing that
##fixed in excel stage legacy code 
##table_input_ordered <- read.csv("4iso - 4iso.csv", sep=",", na.strings = c("","NA"))
##table_input_ordered_df <- as.data.frame(table_input_ordered)
##colnames(table_input_ordered_df) <- colnames(table_input_df)
##View(table_input_ordered_df)


tableValues <- within(table_input, rm("Annotation", "Non.unique.Gene.name","No..isolates","No..sequences","Avg.sequences.per.isolate","Genome.Fragment","Order.within.Fragment","Accessory.Fragment","Accessory.Order.with.Fragment","QC","Min.group.size.nuc","Max.group.size.nuc","Avg.group.size.nuc"))
##View(tableValues)

abscence_presence <- as.matrix(tableValues[,-1])
rownames(abscence_presence) <- tableValues[,1]
abscence_presence[is.na(abscence_presence)] <- 0
abscence_presence[which(abscence_presence!=0)] <- 1

a_p_matrix <- mapply(abscence_presence, FUN=as.numeric)
a_p_matrix <- matrix(data=a_p_matrix, ncol=length(colnames(abscence_presence)), nrow=length(row.names(abscence_presence)))
row.names(a_p_matrix) <- row.names(abscence_presence)
colnames(a_p_matrix) <- colnames(abscence_presence)
##colnames(a_p_matrix) <- c("020MIC", "13121.1", "15108", "1709.12", "2192", "39016", "39177", "57P31PA", "5C1_S24", "679", "968333S", "A5803", "AA2", "AES.1R_2", "AMT0023.30", "AMT0023.34", "AMT0060.1", "AMT0060.2", "AMT0060.3", "AUS058", "AUS066", "AUS088", "AUS089", "AUS23", "AUS52", "C3179", "CHA", "CPHL9433", "DK2", "DUN.001C", "DUN.003B", "DUN.004", "ILPAO1", "IST27", "IST27N", "Jpn1563", "KK1", "LES400", "LES431", "LESB58", "LMG14084", "Mi162_2", "NH57388A","PAO1", "PAK","PAO1", "PAO1_roary", "Pr335", "PAO1", "PA14", "TBCF10839", "U018a", "PA14" )
##View(a_p_matrix)
##View(abscence_presence)

##I want to remove genes that only show up in one strain and genes that show up in all strains
##I have decided to cheat and add a temporary rowsum coloumn and then filter it, and then slice off the extra row
##turns out you cant use $ to subset a matrix, have to use cbind instead 
NewMatrix2 <- cbind(a_p_matrix, newCol = rowSums(a_p_matrix))
View(NewMatrix2)

##I am loosing the will to live, lets make this a dataframe
cheating <-as.data.frame(NewMatrix2)
##View(cheating)
##Removes any genes that only show up once
subsetcheating1 <- cheating[!(cheating$newCol<=1),]
#View(subsetcheating1)
##Removes genes that are in every strain
subsercheating2 <- subsetcheating1[!(subsetcheating1$newCol >=49 ),]
#View(subsercheating2)
##Removes my cheating column
altered_ap <- subsercheating2[,1:50]
#altered_ap <- as.numeric(altered_ap)
#View(altered_ap)
##make this back into a matrix
altered_ap_matrix <- as.matrix(altered_ap)
#View(altered_ap_matrix)
#levels(altered_ap_matrix)

##

heatmap(altered_ap_matrix, Rowv = NA, Colv = NA, labRow=TRUE, cexCol = 0.5,col = c("#FFE986","#FF736E"), par(mar = c(15,4,0.25,4)))


hotpdf<-pdf(hot, width=100, height = 1000)
View(hotpdf)

##This second half I have not touched at all pretty much. Certainly could though!
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
