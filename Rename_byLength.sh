#!/bin/bash

while getopts a:o: option
do
        case "${option}"
        in
        a) Fasta=${OPTARG};;
	o) DIR=${OPTARG};;
        
        esac
done

if [ -z "$Fasta" ] || [ -z "$DIR" ]; then
        echo "
	Usage: /opt/Rename_byLength.sh -a FASTA -o OUTDIR
	
	arguments:
	   -a FASTA
	      FASTA file to be renamed
  	   -o OUTDIR
	      Output directory
	" 1>&2
        exit 1
fi


#Rename scaffold name from reference

FArename=`echo $(basename "$Fasta") | awk '{split($1,s,".fa"); print s[1]}'`
awk '{if(substr($1,1,1)==">")n=$1;else len[n]+=length($1);}END{for(i in len)print i, len[i];}' $Fasta | sort -k2,2nr | awk '{i++; print $1"\t"">scaffold"i"|size"$2}' > $DIR/${FArename}_RenameBylength.tsv
awk 'NR==FNR{scaf[$1]=$2; next} {if(substr($1,1,1)==">")print scaf[$1];else print}' $DIR/${FArename}_RenameBylength.tsv $Fasta > $DIR/${FArename}_rename.fa

