import pandas as pd
import os, glob
import sys
from biopandas.pdb import PandasPdb
import pickle




#Uses argument1 as the path of the folder in all of the .pdb output files
path=str(sys.argv[1])

part=str(path + '/') #string to help format prot accessions for the final dict
my_dict={} #dictionnary in wich i store means
def plddt_mean_calc(filename_): #fonction thats gonna open the pdb select the data that i need, calculate the plddt_mean and then add it to a dictionnary
    plddt_mean=0 #reset plddt_mean

    pdb=PandasPdb().read_pdb(filename)#open the pdb with biopandas
    pdb_unpack=pdb.df['ATOM']#unpack the right biopandas object

    pdb_df=pdb_unpack.loc[:,['atom_name','residue_number','b_factor']] #subset the df with needed values
    pdb_df_ca=pdb_df[pdb_df['atom_name']=='CA'] #subset the df to keep only the alpha carbon lines

    unique_test=pd.Series(pdb_df_ca['residue_number']) # subset for control to make sure the subsetting is alright

    is_unique=unique_test.is_unique #control
    

    #print(filename) #debugstuff
    #print(pdb_df_ca['b_factor'].head()) #debugstuff
    
    if is_unique == True: #makes sure to only execute the fonction if the control conditions are met
        plddt_mean=pdb_df_ca['b_factor'].mean() #mean calculation
        #print(plddt_mean) #debugstuff
        prot_accession=str(filename).partition(part)[2].replace('.pdb','') #format the path to only keep the accession number
        my_dict[prot_accession]=plddt_mean #output dictionnary append of desired values
        return my_dict
    else: #if controls fail warns the user
        print('something went wrong with subsetting')

  


for filename in glob.glob(os.path.join(path,'*.pdb')): #usage of the ffonction for every .pdb file in the specified output folder
    plddt_mean_calc(filename)


with open(str('stuff.pickle'),'wb') as handle: #saves the script output as pickle object
    pickle.dump(my_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    

