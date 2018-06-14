#Katie Noah, June 14, 2018

#set up environment
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
import glob

#bring in all files to be checked for snps
files=glob.glob('*.fasta_aln')

#create blank dataframes to fill in with snps
df_snp=pd.DataFrame()
df_snp_ct=pd.DataFrame()

for file in files:
    print(file)    
    #get sequence information
    seq=SeqIO.parse(file,"fasta")
    first_record=next(seq)
    try:
        second_record=next(seq)
    except:
        continue
    strain=second_record.description.split('|')[0]
    df_snps=pd.DataFrame(columns=['strain','PAO1_pos','PAO1_nt','strain_nt'])
    pos=-1
    rows_list=[]
    inser=''
    delet=''
    temp_pos=0
    start_i=0
    gene=first_record.description.split('|')[-1]

    #look for nucleotides/amino acids that aren't similar between the two strains
    for i in range(len(first_record)):
        if first_record[i]=='-' :
            inser+=second_record[i]
            if first_record[i-1]!='-' or i==0:
                temp_i=i-1
        elif second_record[i]=='-':
            delet+=first_record[i]
            if temp_pos==0:
                temp_pos=pos
            if second_record[i-1]!='-':
                start_i=i-1
            #print(start_i)
            pos+=1
        else:
            if inser!='':
                dict1={}
                inser=second_record[temp_i]+inser
                dict1.update({'strain': strain,'PAO1_pos':pos,'PAO1_nt':first_record[temp_i],'strain_nt':inser})
                rows_list.append(dict1)
                inser=''
            if delet!='':
                dict1={}
                delet=first_record[start_i]+delet
                dict1.update({'strain': strain,'PAO1_pos':temp_pos,'PAO1_nt':delet,'strain_nt':second_record[start_i]})
                rows_list.append(dict1)
                delet=''
                temp_pos=0
            pos+=1
            if first_record[i] != second_record[i]:
                dict1={}
            
                dict1.update({'strain': strain,'PAO1_pos':pos+1,'PAO1_nt':first_record[i],'strain_nt':second_record[i]})
        
                rows_list.append(dict1)
    
    #add to master dataframe        
    try:
        df_snps=df_snps.append(rows_list)
        df_snp_ct=df_snp_ct.set_value(gene,strain,df_snps.shape[0])
        for i, row in df_snps.iterrows():
            if len(row['PAO1_nt'])>1:
                snp=str(row['PAO1_pos']+1)+'_del'
            elif len(row['strain_nt'])>1:
                snp=str(row['PAO1_pos']+1)+'_ins'
            else:
                snp=str(row['PAO1_nt'])+str(row['PAO1_pos'])+str(row['strain_nt'])
            df_snp=df_snp.set_value(gene+'_'+snp,row['strain'],1)
            
    except:
        df_snp_ct=df_snp_ct.set_value(gene,strain,0)
#output files    
df_snp_ct.to_csv('snp_count.csv')
df_snp.to_csv('snps.csv')