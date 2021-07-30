#!/bin/bash
DIR=$PWD
MappedLen=1000
Identity=70

while getopts f:o:a:l:c: option
do
        case "${option}" in
                f) BXCntFile=${OPTARG};;
                o) DIR=${OPTARG};;
                a) Fasta=${OPTARG};;
		l) MappedLen=${OPTARG};;
		c) Identity=${OPTARG};;
        esac
done

if [ -z "$BXCntFile" ] || [ -z "$Fasta" ]; then
        echo "
	usage: /opt/LCB_Scaffolding/Scaffolding.sh
               -f SCAFFOLD_PAIR -a FASTA [-l MIN_MAPLEN] [-c MIN_MAPIDENTITY] -o OUTDIR

        positional arguments:
               -f SCAFFOLD_PAIR
                  Candidate scaffold pairs and shared barcode counts
                  (C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv)
               -a FASTA
                  Draft assembly for LCB scaffolding
               -o OUTDIR
                  Output directory (same as Preprocessing.sh)

        optional arguments:
               -l MIN_MAPLEN
                  Minimum mapped length for contig on each scaffold end [default=1000]
               -c MIN_MAPIDENTITY
                  Mapping identity for contigs on each scaffold end [default=70]
        " 1>&2
	#Usage: $0 [-f: ScaffoldEndPair.tsv] [-a: Fasta file for inter-scaffoling] [-l: Minimum mapped length for contig on each scaffold end] [-c: Mapping identity fot contigs on each scaffold end] [-o: Output directory for assembling step]" 1>&2
        exit 1
fi

SECONDS=0
Rec=$DIR/LCBScaffold.log
printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Start to find contigs crossing scaffold ends..." > $Rec

ConnectInfoFile=$DIR/connect.log
while read -r scaf1 scaf2 BXCnt; do
	s1=`echo $scaf1 | awk '{split($1,a,"_"); print a[1];}'`
        p1=`echo $scaf1 | awk '{print substr($1,length($1)-3,4)}'`
        s2=`echo $scaf2 | awk '{split($1,a,"_"); print a[1];}'`
        p2=`echo $scaf2 | awk '{print substr($1,length($1)-3,4)}'`
	
	wd=$DIR/${s1}${p1}_$s2$p2
	sp1=`echo $scaf1 | awk '{print substr($1,length($1)-3,1)}'`
	sp2=`echo $scaf2 | awk '{print substr($1,length($1)-3,1)}'`


	if [ -f $wd/scaffolds.fasta ]; then
		if [ $sp1 = $sp2 ]; then
			rr="s2_Reverse"
			if [ $sp1 = "T" ]; then
                                sp2="H"
                        else
                                sp2="T"
                        fi
			echo ${s1}${p1}_$s2$p2 s2_Reverse >> $ConnectInfoFile
		else
                        rr="s2_forward"
			echo ${s1}${p1}_$s2$p2 >> $ConnectInfoFile
                fi

		if [ ! -f $wd/Scaf1_aln.sam ]; then
			#Extract scaffold1 and scaffold2 
			python /opt/ExtractScaffoldByName.py $Fasta $s1 > $wd/Scaf1.fa
			python /opt/ExtractScaffoldByName.py $Fasta $s2 > $wd/Scaf2.fa
			if [ $rr = "s2_Reverse" ]; then
				python /opt/scaffold_RC.py $wd/Scaf2.fa > $wd/temp.fa
                        	mv $wd/temp.fa $wd/Scaf2.fa
	        	fi

			#Align contigs to scaffolds
			bwa index $wd/Scaf1.fa
			bwa mem $wd/Scaf1.fa $wd/scaffolds.fasta > $wd/Scaf1_aln.sam

			bwa index $wd/Scaf2.fa
        	        bwa mem $wd/Scaf2.fa $wd/scaffolds.fasta > $wd/Scaf2_aln.sam
		fi

		python /opt/LCB_Scaffolding/CrossContigInfo.py $wd/Scaf1_aln.sam $sp1 $wd/Scaf2_aln.sam $sp2 $MappedLen $Identity >> $ConnectInfoFile
		rm $wd/Scaf1.fa.* $wd/Scaf2.fa.* 		
		echo ============================================= >> $ConnectInfoFile
	fi

done < $BXCntFile

printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Scaffolding with info in $DIR/connect.log..." >> $Rec
python /opt/LCB_Scaffolding/Scaffolding.py $ConnectInfoFile $Fasta

printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "LCB-Scaffolding is done. Collecting final fasta..." >> $Rec
FA=`echo $(basename "${Fasta}") | awk '{split($1,s,".fa"); print s[1]"_LCBScaffold.fa"}'`
if [ -f $DIR/NewScaf.log ];then
        cat $DIR/NewScaf.log | cut -f 1 | xargs -n1 -I {} cat {}.fa > $DIR/$FA
        awk 'NR==FNR{for(i=2;i<=NF;i++)re[">"substr($i,1,length($i)-1)]++; next} {if(substr($1,1,1)==">"){if(re[$1])skip=1;else skip=0;} if(skip==0)print}' $DIR/NewScaf.log $Fasta >> $DIR/$FA
fi
printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "LCB-scaffolding is finished. Concatenated fasta is at $DIR/$FA. Log: $DIR/NewScaf.log." >> $Rec

#Rename scaffolds for next scaffolding
FArename=`echo $(basename "${FA}") | awk '{split($1,s,".fa"); print s[1]}'`
printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Renaming $FA; Renamed Fasta is @$DIR/${FArename}_rename.fa" >> $Rec
awk '{if(substr($1,1,1)==">")n=$1;else len[n]+=length($1);}END{for(i in len)print i, len[i];}' $FA | sort -k2,2nr | awk '{i++; print $1"\t"">scaffold"i"|size"$2}' > $DIR/${FArename}_RenameBylength.tsv
awk 'NR==FNR{scaf[$1]=$2; next} {if(substr($1,1,1)==">")print scaf[$1];else print}' $DIR/${FArename}_RenameBylength.tsv $FA > $DIR/${FArename}_rename.fa
printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Renamed LCB-Scaffolded fasta is at $DIR/${FArename}_rename.fa" >> $Rec

