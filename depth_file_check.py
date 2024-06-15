# This script reads depth files and filters out those files where depth <50
#usage: python check.py -d /path/to/dir/where/Batch60_07012023
#!/usr/bin/env python
import pandas as pd
import os 
import sys
import glob
import argparse
parser = argparse.ArgumentParser("This code compares the dst file results between cryptic and pipeline")
input = parser.add_argument_group(title="Input Files",description = "Enter the required input files")
input.add_argument("-d", "--depthfile",action ='store',dest='dpt', help = "give path of directory where depth files are stored") 
args = parser.parse_args()
Dir_name =args.dpt
 # empty list for reading files
tp=[] 
sf=[] #empty list for storing names
print("reading files.......")

def ffile(Dir_name,fpath):
    rd =[]
    for name in glob.glob(Dir_name+fpath,recursive=True):    # searches for all files in dir and subdir with name ending in _depths.txt
        rd.append(name) 
    return rd
    

def depth():
    rd =ffile(Dir_name,"/*_all_positions_depth.csv")
    if len(rd)<1:
        print("files not found") 
        sys.exit()
    else:                                              #stores path in list rd
        print("Depth files read...............")
        for file in rd: 
            df1 = pd.read_csv(file)                         
            if (df1.V3<50).any() ==True:  
                v=os.path.basename(file).split("_all_positions_depth.csv")[0]     #checks if any value <50 and if true takes basename split at _depths.txt and appends in list
                tp.append(v)
        print("Depth files read \nCreating output textfile......")
    return tp
tp =depth()


def fsummary() :
    sm=ffile(Dir_name,"/*_batch_summary.csv" )
    print(sm)
    if len(sm)<1:
        print("files not found") 
        sys.exit()
    else: 
        print("Reading summary files.....")       
        for file in sm: 
            df1 = pd.read_csv(file)                         
            sf.extend(df1["Sample_name"][df1[" Lineage"]==" unknown"])           #sample name for those whose lineage is unknown                       
            sf.extend(df1["Sample_name"][df1[" MTB_detected"]==" Not detected"]) # sample name for those whose MTB detected is not detected
            print(sf)
    return sf
   
sf =fsummary() 
print("Creating text file")
with open("depth_fail.txt", "a") as fh:
        fh.write("\n".join(tp))
        fh.write("\n")
        fh.write("\n".join(set(sf)))
print("Done")
   
