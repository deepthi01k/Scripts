#!bin/bash

source activate MTBseq
file_name="$(find . -iname "*_R1.fastq.gz" -exec basename {} _R1.fastq.gz ';')"

#move the files with similiar prefix to the folder with the same name that is generated respectively
ls *_R1.fastq.gz | awk -F"_R1.fastq.gz" '{system("mkdir " $1 "/ && mv " $0 " "$1"_R2.fastq.gz " $1 "/" )}'

for ind_sample in $file_name
	do
		cd $ind_sample
		echo "List of paired samples of $ind_sample: "
	   	echo "$ind_sample"_R1.fastq.gz
	     	echo "$ind_sample"_R2.fastq.gz
        	R1="$ind_sample"_R1.fastq.gz
        	R2="$ind_sample"_R2.fastq.gz
	
		eval "MTBseq --step TBfull --threads 12 --lowfreq_vars --minfreq 20 --mincovf 2 --mincovr 2  --resilist /home/alphabox0006/anaconda3/envs/MTBseq/share/mtbseq-1.0.4-2/var/res/old_catalog_resilist.txt --intregions /home/alphabox0006/anaconda3/envs/MTBseq/share/mtbseq-1.0.4-2/var/res/old_catalog_intregions.txt"

    		cd ../
    	done
    	
conda deactivate 
