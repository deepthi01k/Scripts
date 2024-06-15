#this script takes in depth file and annotates positions with genes associated with them
#command : python depth_annotate.py -d <path to depth file> -t <path to extended catalog.txt>
#!/usr/bin/env python
import polars as pl
import argparse
parser = argparse.ArgumentParser("Description")
input = parser.add_argument_group(title="Input Files",description = "Enter the required input files")
input.add_argument("-d", "--depthfile",action ='store',dest='dpt', help = "give path of directory where depth files are stored") 
input.add_argument("-t", "--text",action ='store',dest='ext_cat', help = "give path of extended catalog.txt ")
args = parser.parse_args()
depth_file =args.dpt
extended_cat=args.ext_cat

def read_files(df1,val):
    if val==True:
        depth=pl.read_csv(df1,has_header=False,new_columns=["Genome","Position","depth"],separator="\t")
    else:
        depth=pl.read_csv(df1,separator="\t")
    return depth
depth=read_files(depth_file,True)
catalog=read_files(extended_cat,False)
deep=depth["depth"]
annotations={}

def annotation():
    for rows in depth.rows(named=True):
        position=rows["Position"]
        depth_val=rows["depth"]
        gene_name=""
        for g_row in catalog.rows():
            g_start=g_row[1]
            g_stop=g_row[2]
            g_name=["gene"]
            if position >= g_start and position <= g_stop:
                gene_name = g_row[0]
                break
        annotations[position] = gene_name
    df = pl.DataFrame(data=list(annotations.items()))
    df=df.with_columns(pl.Series(name="Depth",values=deep))
    df.filter(pl.col("column_1").str.lengths()>0).write_csv("Annotated.tsv",separator="\t")
    print("Done")
    return df
annotation()

