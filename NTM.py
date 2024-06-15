# reads the file ntm_log.txt and extracts the microbe name and the values from text file and converts to csv 
# command python NTM.py
import pandas as pd
import re
import warnings
warnings.filterwarnings('ignore')
with open("ntm_log.txt",'r') as fh:   # you can change file location by providing path
     p = (fh.readline()).replace("\n","")
     k = fh.readlines()
     print(len(k))
if len(k) ==0:
     df3 = pd.DataFrame(columns = ["No best hits found"]) 
     df3.to_csv("Ntm_except.csv",index= False)
     print("Else statement executed and Ntm_except Csv Created")
else:
    for count, line in enumerate(k):
        pass

    l = count+1
    #print("length = ",l)
    #print(bool(re.search(p,"this is first line")))
    if ( bool(re.match(p,"this is first line")) and l > 1) == True :
        df2 = pd.DataFrame()
        nlist =[]
        for i in k:
            if ".fasta" in i:
                nlist =[i.replace(".fasta","").replace("\n","")]
            elif bool(re.search(r'\d', i)) == True:
                res = re.findall("\d+\.*\d+",i)
                nlist.extend(res)
       
                # print(nlist)
            if len(nlist)==5: 
                df= pd.DataFrame([nlist],columns = ['Organism','query sample','total genome','query bases','reference bases'])
                nlist=[]
                
                df2 = df2.append(df) #append df2 only after 5 elements in 1 row

        df2["Organism for Best Ntm based on query bases"]= ""
        df2["Organism for Best Ntm based on ref bases"]=""
        maxofquery = df2["query bases"].max()
        maxofref = df2["reference bases"].max()

        def result(base,bestntm,maxval):
            for m in range(len(df2)) :
                 if (df2[m:m+1][base] == df2[base].max()).bool() == True:
                    # print(df2[m:m+1]["Organism"])
                     ntm= list(df2[m:m+1]["Organism"])
                     df2[bestntm][m:m+1] = ntm
                     for line in k:
                        if maxval in line:
                           df2[bestntm][m:m+1]= df2[bestntm][m:m+1]+'\n'+line
                                             
        result("reference bases","Organism for Best Ntm based on ref bases",maxofref)
        result("query bases","Organism for Best Ntm based on query bases",maxofquery)
        df2.to_csv("Manual_Analysis_Results.csv",index= False)
        print("Result csv created")
    
