#!/bin/bash
Thread=16
DIR=$PWD
Distance=30000
ContigLen=1000
ContigCov=2
MappedLen=300
MappedIdentity=80

while getopts n:o:d:t:a:l:c:M:I: option
do
        case "${option}"
        in
        o) DIR=${OPTARG};;
        t) Thread=${OPTARG};;
        n) Num=${OPTARG};;
	a) Fasta=${OPTARG};;
	d) Distance=${OPTARG};;
	l) ContigLen=${OPTARG};;
	c) ContigCov=${OPTARG};;
	M) MappedLen=${OPTARG};;
	I) MappedIdentity=${OPTARG};;
        esac
done

if [ -z "$Fasta" ]; then
        echo "
	usage: /opt/LAB_GapFilling/Fill.sh -a FASTA [-n NUM ] [-t THREADS ] [-l MIN_LEN]
	[-c MIN_COV] [-d MAX_DIS] [-M MIN_MAPLEN] [-I MIN_MAPIDENTITY] -o OUTDIR

        positional arguments:
                   -a FASTA
                      Draft assembly FASTA file to be gap-filled
                   -o OUTDIR
                      Output directory (same as outdir of ProduceBXList.sh)

        optional arguments:
                   -n NUM
                      Top n scaffolds [default=ALL]
                   -t THREADS
                      Number of threads [default=8]
                   -l MIN_LEN
                      Minimum length of assembled contigs [default=1000]
                   -c MIN_COV
                      Minimum coverage of assembled contigs, this property is produced by
		      SPAdes assembler [default=2]
                   -d MAX_DIS
                      Maximum distance between flanking and contigs [default=30000]
                   -M MIN_MAPLEN
                      Minimum mapped length of contigs [default=300]
                   -I MIN_MAPIDENTITY
                      Minimum mapping identity of contigs [default=80]
" 1>&2
	#echo "Usage: $0 [-a: FASTA file for gap-fiiling] [-n: Number of scaffolds to be processed] [-d: Maximum distance between flanking and contigs] [-l: Minimum length of assembled contigs] [-c: Minimum coverage of assembled contigs] [-M: Minimum mapped length of contigs] [-I: Minimum mapping identity of contigs] [-o: Output directory of Assemble.sh] [-t: Number of threads]" 1>&2
        exit 1
fi


FArename=`echo $(basename "$Fasta") | awk '{split($1,s,".fa"); print s[1]}'`
if [ -z "$Num" ]; then
        Num=`wc -l $DIR/${FArename}_RenameBylength.tsv | awk '{print $1}'`
else	
	RestofFA=$DIR/Reference_rename_afterScaf$Num.fa
fi

Filling () {
s=$2
wd=$1/scaf$s
ContigLen=$3
ContigCov=$4
MappedLen=$5
MappedIdentity=$6
Distance=$7

if [ -d $wd ]; then
	#No barcoded reads to fill the gaps
	if [ ! -f $wd/gapRange_BXs.txt ]; then
		echo $wd/scaffold$s.fa >> $1/getFA.txt
		return
	fi 

	#No gap
	g=`wc -l $wd/gap_pos_scaf${s}.txt | awk '{print $1;}'`
	if [ "$g" -eq "0" ]; then
		echo $wd/scaffold$s.fa >> $1/getFA.txt
		return
	fi
	
	gapCnt=`wc -l $wd/gapRange_BXs.txt | awk '{print $1}'`
	python /opt/LAB_GapFilling/getContig_byGap.py $wd $gapCnt $ContigLen $ContigCov > $wd/scaffolds_all.fasta

	#No assembled contigs to fill the gaps
	scafCnt=`ls -l $wd/scaffolds_all.fasta | awk '{print $5}'`
	if [ "$scafCnt" -eq "0" ]; then
		echo $wd/scaffold$s.fa >> $1/getFA.txt
                return
	fi
	
	bwa index $wd/scaffold$s.fa
	mkdir $wd/filled
	bwa mem $wd/scaffold$s.fa $wd/scaffolds_all.fasta > $wd/filled/Scaf${s}_aln.sam
	awk '{r[$1]++; if(r[$1]==1) print;}' $wd/filled/Scaf${s}_aln.sam | cut -f 1,2,4,5,6,16 > $wd/filled/scaf${s}_contig.txt
	python /opt/LAB_GapFilling/contigRange.py $wd/filled/scaf${s}_contig.txt > $wd/filled/scaf${s}_contigRange.txt
	awk -v d=$Distance 'NR==FNR{boundL[$1]=$4-d; boundR[$1]=$5+d; next} {split($1,t,"_"); id=t[7]; if($4>=boundL[id] && $3<=boundR[id])print}' $wd/gapRange_BXs.txt $wd/filled/scaf${s}_contigRange.txt > $wd/filled/scaf${s}_contigRange_byGap.txt

	awk 'NR==FNR{gap[i]=$1;gapEnd[i]=$2;gapC[i]="{";gapPL[i]="{";gapPR[i]="{";i++;next} {s[$1]["complete"]="{"; s[$1]["partial"]="{"; for(i in gap){ \
	if(gapEnd[i]<$4 && gap[i]>$3){ s[$1]["complete"]=s[$1]["complete"]""gap[i]";"; gapC[i]=gapC[i]""$1";";} \
	else if((gap[i]-$4<0?$4-gap[i]:gap[i]-$4)<=50){s[$1]["partial"]=s[$1]["partial"]""gap[i]";"; gapPL[i]=gapPL[i]""$1";";} \
	else if ((gapEnd[i]-$3<0?$3-gapEnd[i]:gapEnd[i]-$3)<=50){s[$1]["partial"]=s[$1]["partial"]""gap[i]";"; gapPR[i]=gapPR[i]""$1";";} } \
	}END{for(j in gap) print gap[j]" "gapEnd[j]" complete: "gapC[j]"} partialL: "gapPL[j]"} partialR: "gapPR[j]"}"}' \
	$wd/gap_pos_scaf$s.txt $wd/filled/scaf${s}_contigRange_byGap.txt | sort -k1,1n > $wd/filled/scaf${s}_coveredContig.txt
		
	python /opt/LAB_GapFilling/filling.py $wd/filled/scaf${s}_coveredContig.txt $wd/filled/scaf${s}_contigRange_byGap.txt $wd/scaffold${s}.fa $wd/scaffolds_all.fasta $MappedLen $MappedIdentity > $wd/filled/filling.log 
	mv $wd/scaffold${s}_filled.fa $wd/filled/Scaf${s}_filled.fa
	echo $wd/filled/Scaf${s}_filled.fa >> $1/getFA.txt
	rm $wd/filled/Scaf${s}_aln.sam $wd/scaffolds_all.fasta $wd/filled/scaf${s}_contig.txt $wd/filled/scaf${s}_contigRange.txt $wd/filled/filling.log $wd/scaffold$s.fa.*

	python /opt/LAB_GapFilling/locateGap.py $wd/filled/Scaf${s}_filled.fa > $wd/filled/gap_pos_scaf${s}_filled.txt

	outRecord=$wd/filled/Scaf${s}_gapFilling.log
	Total_gaps=$(awk '{print $1;}' $wd/gap_pos_scaf$s.txt | wc -l)
	echo 'Total gaps:' $Total_gaps > $outRecord
	Filled_gaps=$(awk '$6=="Filled"{print $1;}' $wd/filled/gapStatus_Record.txt | uniq | wc -l)
	F_per=$(echo $Filled_gaps $Total_gaps | awk '{printf("%.2f", ($1/$2)*100)}')
	UF_per=$(echo $F_per | awk '{printf("%.2f", 100-$1)}')

	echo 'Filled gaps:' $Filled_gaps '('$F_per'%)' >> $outRecord
	echo '---Completely filled gaps:' $(awk '$3=="Complete"{print $1;}' $wd/filled/gapStatus_Record.txt | wc -l) >> $outRecord
	echo '---partially filled gaps:' $(awk '$3!="Complete" && $3!="X" && $3!="Delete"{print $1;}' $wd/filled/gapStatus_Record.txt | uniq | wc -l) >>  $outRecord
	echo '---Delete gaps:' $(awk '$3=="Delete"{print $1;}' $wd/filled/gapStatus_Record.txt | uniq | wc -l) >> $outRecord
	echo ' ' >> $outRecord
	echo 'Unfilled gaps:' $((Total_gaps-Filled_gaps)) '('$UF_per'%)' >> $outRecord
	echo '---Inconsistent gaps:' $(awk '$6=="Inconsistent" || $6=="inDeletion"{print $1;}' $wd/filled/gapStatus_Record.txt | wc -l) >> $outRecord
	echo '---Conflict gaps:' $(awk '$6=="Conflict"{print $1;}' $wd/filled/gapStatus_Record.txt | wc -l) >> $outRecord
	echo '---Unfillale gaps:' $(awk '$6=="noContig"{print $1;}' $wd/filled/gapStatus_Record.txt | wc -l) >> $outRecord
	echo ' ' >> $outRecord
	Total_gapLength=$(awk '{if($4=="")len+=$3; else len+=$4;}END{print len;}' $wd/gap_pos_scaf$s.txt)
	echo "Total Gap Length in Scaffold: " $Total_gapLength " (bp)" >> $outRecord
	Remain_gapLength=$(awk 'BEGIN{len=0};{len+=$3}END{print len;}' $wd/filled/gap_pos_scaf${s}_filled.txt)
	echo "Remaining Gap Length After Gap-Filling: " $Remain_gapLength " (bp)" >> $outRecord
	echo ">>> " $(echo $Total_gapLength $Remain_gapLength | awk '{printf("%.2f", ($1-$2)*100/$1)}') "% of Gaps are filled! <<<" >> $outRecord
#	echo scaf$s Done!!
fi
}
export -f Filling

Rec=$DIR/LABFill_Record.log
SECONDS=0
if [ -f $DIR/getFA.txt ]; then
	rm $DIR/getFA.txt
fi

#Gap filling
printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Filling $Num scaffolds with $Thread threads..." > $Rec
seq 1 $Num | xargs -n 1 -I {} -P $Thread bash -c 'Filling $1 $2 $3 $4 $5 $6 $7' sh $DIR {} $ContigLen $ContigCov $MappedLen $MappedIdentity $Distance

#Collecting fasta
printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Gap filling is done. Start to collect fasta..." >> $Rec
cat $DIR/getFA.txt | xargs -n 1 cat > $DIR/${FArename}_LABFilled.fa
#rm $DIR/getFA.txt
if [ ! -z $Num ] && [ -f $RestofFA ]; then
	cat $RestofFA >> $DIR/${FArename}_LABFilled.fa
#	rm $RestofFA
fi
printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Finished. Gap-filled FASTA is @$DIR/${FArename}_LABFilled.fa" >> $Rec

