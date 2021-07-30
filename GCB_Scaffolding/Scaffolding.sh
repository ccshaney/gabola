#!/bin/bash
DIR=$PWD
MappedLen=1000
Identity=70
Thread=40
Num=2000

while getopts f:o:a:g:t:l:c:n option
do
        case "${option}" in
                f) BXCntFile=${OPTARG};;
                o) DIR=${OPTARG};;
                a) Fasta=${OPTARG};;
		g) contigs=${OPTARG};;
		t) Thread=${OPTARG};;
		l) MappedLen=${OPTARG};;
		c) Identity=${OPTARG};;
                n) Num=${OPTARG};;
        esac
done
if [ -z "$BXCntFile" ] || [ -z "$Fasta" ] || [ -z "$contigs" ]; then
        echo "
	usage: /opt/GCB_Scaffolding/Scaffolding.sh 
               -f SCAFFOLD_PAIR -a FASTA -g G_CONTIGS [-l MIN_MAPLEN]
               [-c MIN_MAPIDENTITY] [-t THREADS] [-n NUM] -o OUTDIR
 
        positional arguments:
                -f SCAFFOLD_PAIR
                   Candidate scaffold pairs and shared barcode counts
                   (C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv)
                -a FASTA 
                   Draft assembly for GCB scaffolding
                -g G_CONTIGS
                   Scaffolds/contigs from other assemblies or long reads for filling in gaps
                -o OUTDIR
                   Output directory (same as Preprocessing.sh)
 
        optional arguments:
                -l MIN_MAPLEN
                   Minimum mapped length for contig on each scaffold end [default=1000]
                -c MIN_MAPIDENTITY
                   Mapping identity for contigs on each scaffold end [default=70]
                -t THREADS
                   Number of tasks/scaffold pairs per run;
                   depending on how many available CPUs [default=40]
                -n NUM
                   Process top n scaffold end pairs [default=2000]

	" 1>&2
	
	#Usage: $0 [-f: ScaffoldEndPair.tsv] [-a: Fasta file for inter-scaffoling] [-s: Contigs for scaffolding][-t: Number of Jobs per run] [-l: Minimum mapped length for contig on each scaffold end] [-c: Mapping identity fot contigs on each scaffold end] [-o: Output directory for assembling step] [-n: process top n pairs]" 1>&2
        exit 1
fi

Scaffolding () {
    	s1=`echo $2 | awk '{split($1,a,"_"); print a[1];}'`
        p1=`echo $2 | awk '{print substr($1,length($1)-3,4)}'`
        s2=`echo $3 | awk '{split($1,a,"_"); print a[1];}'`
        p2=`echo $3 | awk '{print substr($1,length($1)-3,4)}'`

        wd=$1/${s1}${p1}_$s2$p2
	sp1=`echo $2 | awk '{print substr($1,length($1)-3,1)}'`
	sp2=`echo $3 | awk '{print substr($1,length($1)-3,1)}'`
	mkdir $wd
        ConnectLog=$1/${s1}${p1}_${s2}${p2}_subconnect.log 
		if [ $sp1 = $sp2 ]; then
			rr="s2_Reverse"
			if [ $sp1 = "T" ]; then
                                sp2="H"
                        else
                                sp2="T"
                        fi
			echo ${s1}${p1}_$s2$p2 s2_Reverse >> $ConnectLog
		else
                        rr="s2_forward"
			echo ${s1}${p1}_$s2$p2 >> $ConnectLog
                fi

		if [ ! -f $wd/Scaf1_aln.sam ]; then
			#Extract scaffold1 and scaffold2 
			python /opt/ExtractScaffoldByName.py $5 $s1 > $wd/Scaf1.fa
			python /opt/ExtractScaffoldByName.py $5 $s2 > $wd/Scaf2.fa
			if [ $rr = "s2_Reverse" ]; then
				python opt/scaffold_RC.py $wd/Scaf2.fa > $wd/temp.fa
                        	mv $wd/temp.fa $wd/Scaf2.fa
	        	fi

			#Align contigs to scaffolds
			minimap2 -d $wd/Scaf1.fa.fai $wd/Scaf1.fa
			minimap2 -ax asm5 -o $wd/Scaf1_aln.sam $wd/Scaf1.fa $6 2&>$wd/Scaf2_aln.log

			minimap2 -d $wd/Scaf2.fa.fai $wd/Scaf2.fa
                 	minimap2 -ax asm5 -o $wd/Scaf2_aln.sam $wd/Scaf2.fa $6 2&> $wd/Scaf2_aln.log
		fi

		python /opt/GCB_Scaffolding/GCB_CrossContigInfo.py $wd/Scaf1_aln.sam $sp1 $wd/Scaf2_aln.sam $sp2 $7 $8 >> $ConnectLog
		rm $wd/Scaf1.fa.* $wd/Scaf2.fa.* 		
		echo ============================================= >> $ConnectLog
	
}
export -f Scaffolding

SECONDS=0
Rec=$DIR/GCB_Scaffold.log

printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Renaming g-scaffolds..." >> $Rec
CONrename=`echo $(basename "$contigs") | awk '{split($1,s,".fa"); print s[1]}'`
awk '{if(substr($1,1,1)==">")n=$1;else len[n]+=length($1);}END{for(i in len)print i, len[i];}' $contigs | sort -k2,2nr | awk '{i++; print $1"\t"">scaffold"i"|size"$2}' > $DIR/${CONrename}_RenameBylength.tsv
awk -v dir=$DIR 'NR==FNR{scaf[$1]=$2; next} {if(substr($1,1,1)==">"){print scaf[$1]; split(scaf[$1],t,"|"); sidx=substr(t[1],2); p=scaf[$1];}else{print; p=$1} print p > dir"/"sidx".fa";}' $DIR/${CONrename}_RenameBylength.tsv $contigs > $DIR/${CONrename}_rename.fasta

printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Start to find g-scaffolds crossing scaffold ends..." > $Rec
head -n $Num $BXCntFile | xargs -n 1 -I {} -P $Thread bash -c 'Scaffolding $1 $2 $3 $4 $5 $6 $7 $8' sh $DIR {} $Fasta $DIR/${CONrename}_rename.fasta $MappedLen $Identity 
ConnectInfoFile=$DIR/connect.log
ls $DIR/*_subconnect.log | xargs -n 1 cat > $ConnectInfoFile
rm $DIR/*_subconnect.log

printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "GCB-Scaffolding with info in $DIR/connect.log..." >> $Rec
python /opt/GCB_Scaffolding/GCB_Scaffolding.py $ConnectInfoFile $Fasta $contigs

printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "GCB-Scaffolding is done. Collecting final fasta..." >> $Rec
FA=`echo $(basename "${Fasta}") | awk '{split($1,s,".fa"); print s[1]"_GCBScaffold.fa"}'`
if [ -f $DIR/NewScaf.log ];then
        cat $DIR/NewScaf.log | cut -f 1 | xargs -n 1 -I {} cat {}.fa > $DIR/$FA
        awk 'NR==FNR{for(i=2;i<=NF;i++)re[">"substr($i,1,length($i)-1)]++; next} {if(substr($1,1,1)==">"){if(re[$1])skip=1;else skip=0;} if(skip==0)print}' $DIR/NewScaf.log $Fasta >> $DIR/$FA
fi
printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "GCB-scaffolding is finished. Concatenated fasta is at $DIR/$FA. Log: $DIR/NewScaf.log." >> $Rec

#Rename scaffolds for next scaffolding
FArename=`echo $(basename "${FA}") | awk '{split($1,s,".fa"); print s[1]}'`
printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Renaming $FA; Renamed Fasta is @$DIR/${FArename}_rename.fa" >> $Rec
awk '{if(substr($1,1,1)==">")n=$1;else len[n]+=length($1);}END{for(i in len)print i, len[i];}' $FA | sort -k2,2nr | awk '{i++; print $1"\t"">scaffold"i"|size"$2}' > $DIR/${FArename}_RenameBylength.tsv
awk 'NR==FNR{scaf[$1]=$2; next} {if(substr($1,1,1)==">")print scaf[$1];else print}' $DIR/${FArename}_RenameBylength.tsv $FA > $DIR/${FArename}_rename.fa
printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Renamed GCB-Scaffolded fasta is at $DIR/${FArename}_rename.fa" >> $Rec


