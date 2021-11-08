library(DESeq2)

#Read in count information. 
countData = read.table("HamsterCounts.tsv",header=TRUE,row.names=1,sep = '\t')
##Read in expermental design
metadata = read.table("Design",header=TRUE,row.names=1,sep = '\t')




#Get list of tissues 
Tissue=list(metadata$TISSUE)
Tissue= sapply(Tissue, unique)

for (tissue in Tissue) {


    #Select Tissue
    colData=metadata[metadata$TISSUE==tissue,]
    
    
    #Get sample names 
    samples=row.names(colData)
    #Extract samples from counts
    counts=countData[,samples]
    #Round to nearest int
    counts=round(counts,0)
    
    
    
    #Should return TRUE
    #all(rownames(colData) == colnames(counts))
    
    ##Make DEseq2 object 
    dds = DESeqDataSetFromMatrix(countData = counts,colData = colData,design = ~ INFECTION)   
    dds = DESeq(dds)
    #Contrast SARS-CoV-2 vs Mock
    result = results(dds, contrast=c("INFECTION","SARS-CoV-2","Mock"))
    ## Remove rows with NA
    result = result[complete.cases(result),]
    #Put GeneID as column 
    result = cbind(GeneID = rownames(result), result)
    write.table(result, paste("SARS-CoV-2_vs_Mock_", tissue, ".txt",sep=""),sep = '\t',row.names = FALSE)

}
