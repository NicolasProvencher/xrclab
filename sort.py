import numpy as np
import pandas as pd
head=['Gene stable ID','Gene Synonym','Gene name'] #first document head
headdf=['Gene stable ID','Gene Synonym','Gene name','HGNC ID','HGNC symbol','NCBI gene (formerly Entrezgene) accession'] #all the document head added
df=pd.DataFrame(columns=headdf) #creates empty dataframe in wich im gonna add line


with open('ens_gene_synonymes.txt','r') as gene_1: #opens first text doc
    next(gene_1) #skip first line
    for lines in gene_1: #for each entry
        
        ls=lines.strip('\n').split('\t') #format the entry to get a list
        if ls[0] in df[head[0]].values: #if the gene stable id is already in the data frame, 
            x=df[df[head[0]]==ls[0]].index #store retain the index of that gene stable id
            for i,j in enumerate(ls[1:],1): #for all element in the list except the gene stable id
                if j in df.loc[x[0]][head[i]]: #if the element is already there skip
                    pass
                else: #otherwise adds the new element after a ,

                    df.loc[x[0]][head[i]]=df.loc[x[0]][head[i]]+','+j

                    
        else: #thats if the gene stable id isnt present in the data frame, we add that line to the data frame
            df1=pd.DataFrame([ls],columns=head) #create a line to concat
            df=pd.concat([df,df1],ignore_index=True) #concat the line while ignoring index so we can generate a new index

df = df.replace(np.nan,'') #replaces all numpy nun with empty spot so we can eiitherate into later (so one of the if works)
print("1e document done") #anoucne the 1st step done



with open('mart_export_hgnc_genes.txt','r') as gene_2: #opens 2nd document
    head1 = gene_2.readline() #add the first one to a header object
    head1=head1.strip('\n').split('\t') #format the header into workable thing
    
    for lines in gene_2: #pass thorught each line mostly same thing as last document except 1 added if
        
        ls=lines.strip('\n').split('\t')
        if ls[0] in df[head1[0]].values:
            x=df[df[head1[0]]==ls[0]].index

            for i,j in enumerate(ls[1:],1):

                if j in df.loc[int(x[0])][head1[i]]:

                    pass
                elif df.loc[int(x[0])][head1[i]]=="": #new if so that f the space is empty a value is added without a , in front
                    df.loc[int(x[0])][head1[i]]=j
                else:

                    df.loc[int(x[0])][head1[i]]=df.loc[int(x[0])][head1[i]]+','+j
        else:
            df1=pd.DataFrame([ls],columns=head1)
            df=pd.concat([df,df1],ignore_index=True)
print("2e document done")

df = df.replace(np.nan,'')
with open('mart_export_ncbi_genes.txt','r') as gene_3: #same thing as other one
    head3 = gene_3.readline()
    head3=head3.strip('\n').split('\t')
    print(head3)
    
    for lines in gene_3:
        ls=lines.strip('\n').split('\t')
        if ls[0] in df[head3[0]].values:
            x=df[df[head3[0]]==ls[0]].index

            for i,j in enumerate(ls[1:],1):

                if j in df.loc[int(x[0])][head3[i]]:

                    pass
                elif df.loc[int(x[0])][head3[i]]=="":
                    df.loc[int(x[0])][head3[i]]=j
                else:

                    df.loc[int(x[0])][head3[i]]=df.loc[int(x[0])][head3[i]]+','+j
        else:
            df1=pd.DataFrame([ls],columns=head3)
            df=pd.concat([df,df1],ignore_index=True)


df = df.replace(np.nan,'')
print("3e document done")

with open('hgnc_biomart.txt','r') as gene_4: #same thing as other one but with added step to add a number to row with missing ensmb
    head4 = gene_4.readline()
    head4=head4.strip('\n').split('\t')
    head4[-1]=head3[0] #here we use -1 cause the gene stable id is the last column in the document and this line is used to rename it
    print(head4)
    
    for flag,lines in enumerate(gene_4):
        ls=lines.strip('\n').split('\t')
        if ls[-1]=='':
            ls[-1]=flag

        if ls[-1] in df[head4[-1]].values:
            x=df[df[head4[-1]]==ls[-1]].index

            for i,j in enumerate(ls[:-1],0):

                if j in df.loc[int(x[0])][head4[i]]:

                    pass
                elif df.loc[int(x[0])][head4[i]]=="":
                    df.loc[int(x[0])][head4[i]]=j
                else:

                    df.loc[int(x[0])][head4[i]]=df.loc[int(x[0])][head4[i]]+','+j
        else:
            df1=pd.DataFrame([ls],columns=head4)
            df=pd.concat([df,df1],ignore_index=True)








df.to_csv("testingfin.tsv", sep='\t')
