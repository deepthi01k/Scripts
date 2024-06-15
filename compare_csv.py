#this script compares two dataframes of equal shape i.e number of rows and columns are same and have same column names and highlights the differences
#requirement :pip install openpyxl
#!/usr/env python
import pandas as pd
import numpy as np
import openpyxl
import os
import argparse
parser = argparse.ArgumentParser("This code compares 2 csv files and returns unmatched rows in a csv")
input = parser.add_argument_group(title="Input Files",description = "Enter the required input files")
input.add_argument("-f1", "--file1",action ='store',type ='string',dest='csv1', help = "the original csv")
input.add_argument("-f2", "--file2",action ='store',type ='string',dest='csv2', help = "the other csv")
args = parser.parse_args()
name1 = args.csv1
name2 = args.csv2
df1 =pd.read_csv(name1)
df2 =pd.read_csv(name2)
def compare_excel(): 
     df_all = pd.concat([df1, df2], axis='columns', keys=['Set1', 'Set2'])
     df_final = df_all.swaplevel(axis='columns')[df1.columns[1:]]
     def highlight_diff(data, color='lightcoral'):
          attr = 'background-color: {}'.format(color)
          other = data.xs('Set1', axis='columns', level=-1)
          return pd.DataFrame(np.where(data.ne(other, level=0), attr, ''),
                        index=data.index, columns=data.columns)

     compared_csv=df_final.style.apply(highlight_diff, axis=None).to_excel("output.xlsx") #output is an excel sheet
     return compared_csv
compare_excel()
