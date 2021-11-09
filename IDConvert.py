import pandas as pd
import glob

#Read in list of GeneIDs and corresponding GeneNames.
GeneID_GeneName = pd.read_csv("GeneID_GeneName.tsv",sep='\t')

#DGE files
myFilesPaths = glob.glob(r'*.txt')

for file in myFilesPaths:
    #Read in DGE file
    dgeFile = pd.read_csv(file,sep='\t')   
    #Merge on GeneID
    df = pd.merge(GeneID_GeneName,dgeFile,on=['GeneID'],how="outer")
    #Remove Na rows 
    df = df[df['padj'].notna()]
    #write 
    
    
    df.to_csv(file.replace(".txt", ".tsv"),mode='w', header=True,index=False,sep='\t')
