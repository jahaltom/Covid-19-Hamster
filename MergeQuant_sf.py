import pandas as pd
import glob


#List quant.sf files
myFilesPaths = glob.glob(r'*.sf')
#Make empty df
counts = pd.DataFrame()
#Get Gene IDs
Geneids = pd.read_csv(myFilesPaths[1],sep='\t')
Geneids=Geneids[['Name']]


for file in myFilesPaths:
    quant = pd.read_csv(file,sep='\t')
    #Get counts
    quant=quant[['NumReads']]
    #Rename counts column sample name 
    quant.columns=[file.replace(".quant.genes.sf", "")]  
    #Append
    counts = pd.concat([counts, quant], axis=1)
    
#Add gene ID column    
counts = pd.concat([Geneids, counts], axis=1)   

counts.to_csv("HamsterCounts.tsv",mode='w', header=True,index=False,sep="\t")
