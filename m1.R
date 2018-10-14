library(GEOquery)

###############################################################
## Acquire the Microarray Data for Strains AKR/J and C3H/HeJ ##
###############################################################

## Acquire the microarray platform info.
gpl = getGEO("GPL8759")
gpl = Table(gpl)
colnames(gpl)
gpl[1:10, c(1,11)]

## Acquire the GSM data of control cells from the AKR/J strain.
## This data has been RMA-normalized (log2 scale).
control_akr = getGEO("GSM948046")
control_akr = Table(control_akr)
dim(control_akr)
control_akr[1:10,]

## Acquire the GSM data of LPS-treated cells from the AKR/J strain.
## This data has been RMA-normalized (log2 scale).
lps_akr = getGEO("GSM948224")
lps_akr = Table(lps_akr)
dim(lps_akr)
lps_akr[1:10,]

## Acquire the GSM data of control cells from the C3H/HeJ strain.
## This data has been RMA-normalized (log2 scale).
control_c3h = getGEO("GSM948017")
control_c3h = Table(control_c3h)
dim(control_c3h)
control_c3h[1:10,]

## Acquire the GSM data of LPS-treated cells from the C3H/HeJ strain.
## This data has been RMA-normalized (log2 scale).
lps_c3h = getGEO("GSM948194")
lps_c3h = Table(lps_c3h)
dim(lps_c3h)
lps_c3h[1:10,]

#################################################
## Add Gene Symbols to the Microarray Datasets ##
#################################################

control_akr$GENE = paste(gpl[["Gene Symbol"]][match(control_akr$ID_REF, gpl$ID)])
control_akr[1:10,]
lps_akr$GENE = paste(gpl[["Gene Symbol"]][match(lps_akr$ID_REF, gpl$ID)])
lps_akr[1:10,]
control_c3h$GENE = paste(gpl[["Gene Symbol"]][match(control_c3h$ID_REF, gpl$ID)])
control_c3h[1:10,]
lps_c3h$GENE = paste(gpl[["Gene Symbol"]][match(lps_c3h$ID_REF, gpl$ID)])
lps_c3h[1:10,]

###################
## Data Clean-Up ##
###################

## Some probes have multiple names representing multiple genes, such as the 10th row here:
control_akr[1:10,]
## I prefer to reduce these names to their more recognizable singular gene symbols,
## by removing everything before and including the "///" characters.
control_akr$GENE= gsub("^.*?///", "", control_akr$GENE)
lps_akr$GENE= gsub("^.*?///", "", lps_akr$GENE)
control_c3h$GENE= gsub("^.*?///", "", control_c3h$GENE)
lps_c3h$GENE= gsub("^.*?///", "", lps_c3h$GENE)

## Remove the rows with blank values
control_akr[control_akr == ""] = NA
control_akr = control_akr[-which(is.na(control_akr$GENE)),]
which(is.na(control_akr))

lps_akr[lps_akr == ""] = NA
lps_akr = lps_akr[-which(is.na(lps_akr$GENE)),]
which(is.na(lps_akr))

control_c3h[control_c3h == ""] = NA
control_c3h = control_c3h[-which(is.na(control_c3h$GENE)),]
which(is.na(control_c3h))

lps_c3h[lps_c3h == ""] = NA
lps_c3h = lps_c3h[-which(is.na(lps_c3h$GENE)),]
which(is.na(lps_c3h))

######################################
## Differential Expression Analysis ##
######################################

## Subtract control values from LPS values
akr_lps_over_control = within(merge(lps_akr, control_akr, by = "ID_REF"), {VALUE = lps_akr$VALUE-control_akr$VALUE})
akr_lps_over_control = akr_lps_over_control[, -c(2:4)]
colnames(akr_lps_over_control) = c("ID_REF", "GENE", "VALUE")
akr_lps_over_control[1:5,]

c3h_lps_over_control = within(merge(lps_c3h, control_c3h, by = "ID_REF"), {VALUE = lps_c3h$VALUE-control_c3h$VALUE})
c3h_lps_over_control = c3h_lps_over_control[, -c(2:4)]
colnames(c3h_lps_over_control) = c("ID_REF", "GENE", "VALUE")
c3h_lps_over_control[1:5,]

## Remove any duplicate probes (I prefer to take only the top scoring probe of any duplicates).
akr_lps_over_control = akr_lps_over_control[!duplicated(akr_lps_over_control$GENE), ]
dim(akr_lps_over_control)
c3h_lps_over_control = c3h_lps_over_control[!duplicated(c3h_lps_over_control$GENE), ]
dim(c3h_lps_over_control)

###################################################
## Determine Expression of M1 Macrophage Markers ##
###################################################

## Before proceeding, acquire the list of M1 marker genes from the 2010 paper by Kadl et al. (https://www.ncbi.nlm.nih.gov/pubmed/20651288).
## Here, I have named my M1 marker file as "M1 Table.csv".
m1 = read.csv("M1 Table.csv", header = TRUE, stringsAsFactors = FALSE)

## Subset the microarray according to the M1 marker list.
m1_akr = akr_lps_over_control[akr_lps_over_control$GENE %in% m1$Full_M1_Gene_Name,]
dim(m1_akr)
m1_c3h = c3h_lps_over_control[c3h_lps_over_control$GENE %in% m1$Full_M1_Gene_Name,]
dim(m1_c3h)

## Combine the m1_akr and m1_c3h datasets together.
## Then convert the dataset to a matrix. This will be used to create the heat map.
m1_akr_c3h = merge(m1_akr, m1_c3h, by.x = "ID_REF", by.y = "ID_REF")
rownames(m1_akr_c3h) = m1_akr_c3h$GENE.x
m1_matrix = as.matrix(m1_akr_c3h[,-c(1,2,4)])
colnames(m1_matrix) = c("AKR/J", "C3H/HeJ")
head(m1_matrix)

## You can find the max and min values of the matrix.
max(m1_matrix)
min(m1_matrix)

## You can calculate the percentage of M1 marker genes that are upregulated.
## For 1.5 fold change, calculate 2^x = 1.5, which uses the R code: log(1.5, base = 2)
## For the AKR/J strain:
sum(m1_akr$VALUE >= 0.585)/length(m1_akr$VALUE)*100
## For the C3H/HeJ strain:
sum(m1_c3h$VALUE >= 0.585)/length(m1_c3h$VALUE)*100

## You can calculate the percentage of M1 marker genes that are downregulated.
## For less than 1 fold change, calculate 2^x = 1, which uses the R code: log(1, base = 2)
## For the AKR/J strain:
sum(m1_akr$VALUE < 0)/length(m1_akr$VALUE)*100
## For the C3H/HeJ strain:
sum(m1_c3h$VALUE < 0)/length(m1_c3h$VALUE)*100

## Create a heat map to visualize the M1 gene expression for both strains.
## This heat map will have the gene expressions in descending order of the AKR/J strain.
library(gplots)
order(m1_matrix[,1])
m1_matrix_ordered = m1_matrix[order(-m1_matrix[,1]),]
head(m1_matrix_ordered)
## Then create the heat map.
colors = c(seq(-7.3,-0.585, length = 250), seq(-0.584,0.584, length = 100), seq(0.585,7.3, length = 250))
my_palette = colorRampPalette(c("red", "black", "green"))(n = 599)
heatmap.2(m1_matrix_ordered, Rowv = FALSE, Colv = FALSE,
          dendrogram = "none", symm = FALSE, scale = "none",
          breaks = colors, symbreaks = TRUE, col = my_palette, 
          trace = "none", margins = c(8.2,5.5), cexRow = 0.1, cexCol = 1.4,
          offsetCol = -0.1, density.info = "none", symkey = FALSE,
          key.title = NULL,
          key.xlab = expression(paste(~log[2]~(ratio))))
