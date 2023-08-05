BiocManager::install("DESeq2")#one time
library(DESeq2)#each time you open your R

#loading the gene expression data 
b_data = as.matrix(read.csv("bladder_counts.csv", row.names=1))
#loading the phenotype data 
b_pheno = read.csv("bladder_sample_sheet.csv", row.names=1)

table(b_pheno$Sample.Type)

#explore the data, dim function give you the dimension of your data; how 
#many columns(samples) do you have and how many row(genes) do you have
dim(b_data)

#explore the data distribution using the histogram plot
hist(b_data, col = "orange", main="Histogram", breaks = 100)

#scaling the data using log2 transformation to better visulization
# we use (+1) to avoid the infinity character when we log zero valus 
hist(log2(b_data+1), col = "orange", main="Histogram")

#It is absolutely critical that the columns of the "b_data" and the rows of 
#the "b_pheno" (information about samples) are in the same order.DESeq2 will
#not make guesses as to which column of the count matrix belongs to which row
#of the column data, these must be provided to DESeq2 already in consistent order
b_pheno=b_pheno[colnames(b_data),]

#The deseq2 package require the count data values to be integers 
#save the gene names in a variable
genes=row.names(b_data)
#convert the data values to integers
data=apply(b_data,2,as.integer)
#view the data
head(b_data)
#rename the rows of the data
row.names(b_data)=genes
#view the data

###### DO the differential EXP analysis using DeSeq2

#specify how many conditions do you want to compare according to 
#the phenotypic table
condition1="Primary Tumor" 
condition2="Solid Tissue Normal"

#creat a deseq dataset object
dds_bladder= DESeqDataSetFromMatrix( countData = b_data , colData = b_pheno, design = ~ Sample.Type)
#run the deseq2 worflow
dds_b.run = DESeq(dds_bladder)
#specifying teh contrast (to make a res object based on two specific conditions)
res=results(dds_b.run, contrast = c("Sample.Type",condition1 ,condition2))

# remove nulls
res=as.data.frame(res[complete.cases(res), ])

#chose the statstical significant differentaily expressed genes (DEGs) based
#on the p adjusted value less than 0.05 and biological significance  based
#on the fold change more than 2
deseq.deg=res[res$padj < 0.05 & abs(res$log2FoldChange)>2,]
dim(deseq.deg)

#export the Degs into your current folder for further analysthis
write.csv(as.matrix(deseq.deg),file="baldder_deseq.deg.csv", quote=F,row.names=T)

#using nonprametric test as wilcoxon test if the data isn't normally distributed
w_b=row_wilcoxon_twosample(b_data[,1:19],b_data[,20:409])

#correct the wilcoxon test p value using the FDR method
p.adj_w_b=p.adjust(w_b$pvalue, method = "fdr")

#save the result in a dataframe contain the fold change and p adjusted value
result_w_b=as.data.frame(cbind(fold , p.adj_w_b))
dim(result_w_b)

#chose the statistical significant differentaily expressed genes (DEGs) based
#on the p adjusted value less than 0.05 and biological significance  based
#on the fold change more than 2
res.deg_w_bladder=result[result_w_b$p.adj < 0.05 & abs(result_w_b$fold)>2,]
dim(res.deg_w_bladder)

#export the DEGs into your current folder for further analysis
write.csv(as.matrix(res.deg_w_bladder),file="res_baldder_wTest.degs.csv", quote=F,row.names=T)




