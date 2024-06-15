#This script creates indel catalog using who catalog and mycobrowser tsv for H37Rv
# USage : python Indelcat.py  N.B keep the Gen , Mut and Myco in same directory from where u are running the script
#!/usr/bin/env python3
import pandas as pd
import re
import warnings
warnings.filterwarnings("ignore")

print("reading")
# reading the input files
Gen = pd.read_table("Genome_indices.tsv", low_memory=False)
Mut = pd.read_table("Mutation_catalogue.tsv", low_memory=False)
Myco =pd.read_table("Mycobacterium_tuberculosis_H37Rv.tsv", low_memory=False)

def variants():
    df5 =sort_mutation()
    dels =df5[df5["final_annotation.TentativeHGVSNucleotidicAnnotation"].str.contains(pat = r"del*(?!.*ins[AGCT])", regex = True, na=False)].reset_index(drop=True) #negative lookahead for deletion only
    ins =df5[df5["final_annotation.TentativeHGVSNucleotidicAnnotation"].str.contains(pat = r"(?<![AGCT])(ins*)", regex = True, na=False)].reset_index(drop=True) #negative lookbehind for insertion only
    dups_ins=df5[df5["final_annotation.TentativeHGVSNucleotidicAnnotation"].str.contains(pat = r"dup*(?=.*ins[AGCT])", regex = True, na=False)].reset_index(drop=True)#positive lookahead for ins+dups
    indel=df5[df5["final_annotation.TentativeHGVSNucleotidicAnnotation"].str.contains(pat = r"del*(?=.*ins[AGCT])", regex = True, na=False)].reset_index(drop=True)#positive lookahead for indels
    dups =df5[df5["final_annotation.TentativeHGVSNucleotidicAnnotation"].str.contains(pat = r"dup*(?!.*ins[AGCT])", regex = True, na=False)] .reset_index(drop=True)#negative lookahaead for only dups
    insertion=ins[~ins["final_annotation.TentativeHGVSNucleotidicAnnotation"].isin(dups_ins["final_annotation.TentativeHGVSNucleotidicAnnotation"])].reset_index(drop=True)
    return dels,insertion,dups

def indel_type(indel):
     nopat =indel[~indel["final_annotation.TentativeHGVSNucleotidicAnnotation"].str.contains(pat = r",", regex = True, na=False)].reset_index(drop=True)# has no commas
     left=nopat["final_annotation.TentativeHGVSNucleotidicAnnotation"].str.extract(r'(-?\d+)',expand=False)
     right=nopat["final_annotation.TentativeHGVSNucleotidicAnnotation"].str.extract(r'(\_-?\d+)',expand=False).str.replace("_","")
     nopat["Variant position gene start"] =left
     nopat["Variant position gene stop"] =right
     nopat["Variant position gene stop"].fillna(inplace=True, value=nopat["Variant position gene start"])
     nopat["temp"]=nopat["final_annotation.TentativeHGVSNucleotidicAnnotation"].str.extract(r'([AGTC]+)',expand =False)
     nopat["Number"] =nopat["temp"].str.len()
     nopat.drop(["temp"],axis=1,inplace=True)
     if (nopat["final_annotation.TentativeHGVSNucleotidicAnnotation"].str.contains(r"(del)",regex=True)).all()==True:
          nopat["Variant position genome stop"]=nopat["Number"]+nopat["genome_index"].astype(int)
     else:
          nopat["Variant position genome stop"]="-"
     return nopat

def sort_Myco(Myco):
    #   Renaming the columns in order to help in merging
    Myco.rename(columns={"Name":"gene_name"},inplace = True)


    # selecting required columns from Mycobacterium_tuberculosis_H37Rv.tsv
    Myco_select = pd.DataFrame()
    Myco_select= Myco.filter(items =["Start","Stop","gene_name","Strand"])

    # calculating gene length 
    Myco_select['Gene length'] = Myco_select.apply(lambda x: x['Stop'] - x['Start']+1, axis=1) 
    print("Successfully edited Myco_select")
    return Myco_select

def sort_mutation():
    # Filtering genomic indices and merging with grouped mutation catalogue
    Mut.rename(columns={"variant (common_name)":"variant"},inplace = True)
    # sort mutation_catalogue and merge with genomic indices
    df2 = Mut[(Mut["FINAL CONFIDENCE GRADING"] == "1) Assoc w R") |(Mut["FINAL CONFIDENCE GRADING"] =="2) Assoc w R - Interim" )]
    df3 = df2.groupby("variant")["drug"].apply(lambda tags: ' , '.join(tags)).reset_index(name="Antibiotic")
    df5 = df3.merge(Gen, how = "left",on = ["variant"])
    # df5 = df4[df4["final_annotation.TentativeHGVSNucleotidicAnnotation"].str.contains(pat = "[>]", regex = True, na=False)]
    print("Successfully merged df1 and mut")
    
    return df5



print("next function")    
def sort_Gen():
    dels,insertion,dups =variants()     
    inss=indel_type(insertion)
    dele =indel_type(dels)
    dupl =indel_type(dups)
    Deli=pd.concat([inss,dele,dupl],ignore_index=True)
    Myco_select = sort_Myco(Myco) # calling value from function
    
    Deli["alt_aa"]=Deli["alt_aa"].replace(['Ter'],value="STOP") 

    #  creating col codon change
    Deli['Codon change']=['/'.join(i) for i in zip(Deli['ref_nt'].astype(str),Deli['alt_nt'].astype(str))]

    # creating col AA change
    Deli["AA change"] =Deli["ref_aa"]+Deli["codon_number"]+Deli["alt_aa"] 

    df6 =pd.merge(Deli,Myco_select,"left","gene_name")
    # extracting Wt nt and Var nt from HGVS column 
    
    df6["codon_number"].fillna("-",inplace=True)
    # df6.to_csv("genomic.csv")
    filter =pd.DataFrame()
    filter=df6.loc[:,["genome_index","Variant position genome stop",'ref_nt',"codon_number",'alt_nt',"gene_name"]].apply(lambda row: '/'.join(row.values.astype(str)), axis=1)
    df6.insert(0,"refcol",filter)
    # print(df6["final_annotation.TentativeHGVSNucleotidicAnnotation"])
    # Filling in values for variant position gene start and stop
    Vary =list()
    for i in range(len(df6)):
       if bool(re.findall(r"(del)",df6["final_annotation.TentativeHGVSNucleotidicAnnotation"][i]))==True:
          Vary.append("Del")
       elif bool(re.findall(r"(ins)",df6["final_annotation.TentativeHGVSNucleotidicAnnotation"][i]))==True: 
           Vary.append("Ins")
          
       elif bool(re.findall(r"(dup)",df6["final_annotation.TentativeHGVSNucleotidicAnnotation"][i]))==True: 
           Vary.append("dup")
            
        
        
        
    # Filing in values for region
    Type_new = list()
    for i in range(len(df6)):
        if bool(re.findall("(^c\.-)",df6["final_annotation.TentativeHGVSNucleotidicAnnotation"][i]))==True:
        #print(i)
            Type_new.append("5' Upstream")
 
        elif bool(re.findall("(^n\.)",df6["final_annotation.TentativeHGVSNucleotidicAnnotation"][i]))==True:
            Type_new.append("non coding")
        
        else:
            Type_new.append("coding")

     
    df7=df6[['refcol','genome_index','Variant position genome stop','Number','ref_nt','alt_nt','gene_locus','gene_name','Start','Stop','Gene length','Strand','ref_aa','codon_number', 'alt_aa','AA change','Codon change',"Variant position gene start","Variant position gene stop",'Antibiotic']]
    # renaming the drugs according to the Old catalog
    vals_to_replace ={"AMK":"amikacin (AMK)","BDQ":"bedaquiline (BDQ)","CAP":"capreomycin (CAP)","DLM":"delamanid (DLM)","EMB":"ethambutol (EMB)","ETH":"ethionamide (ETH)","INH":"isoniazid(INH)","KAN":"kanamycin(KAN)","LFX":"fluoroquinolones (FQ)","LZD": "linezolid(LZD)","PZA":"pyrazinamide(PZA)","RIF":"rifampicin (RIF)","STM":"streptomycin (STM)","AMI , KAN":"amikacin (AMK) kanamycin (KAN)","LEV , MXF" :"fluoroquinolones (FQ)","AMI , CAP , KAN" :"amikacin (AMK) kanamycin (KAN) capreomycin (CPR)","MFX":"fluoroquinolones (FQ)"}
    df7['Antibiotic'] = df7['Antibiotic'].map(vals_to_replace)
    #    renaming columns
    df7.rename(columns = {'genome_index':'Variant position genome start','gene_locus':'Gene ID','gene_name':'Gene Name','Start':'Gene start','Stop':'Gene stop','Strand':'Dir.','ref_aa':'WT AA','codon_number':'Codon nr.','alt_aa':'Var. AA','ref_nt':'WT','alt_nt':'Var. base'},inplace =True)
    # inserting columns 
    df7.insert(3,"Var. type",Vary)
    # df7.insert(4,"Number",1)
    df7.insert(7, "Region", Type_new)
    df7.insert(22,"Reference PMID","-")
    df7.insert(23,"High Confidence SNP","-")
    df7.insert(16,"Codon nr. E. coli","-")
    
   
    
    df7["refcol"] =df7["refcol"]+'/'+df7["Antibiotic"]
    df7.to_csv("TBCatalog.csv", index =False)
sort_Gen()
print("Done!")




