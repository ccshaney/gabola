#!/bin/bash
set -e

Thread=16
DIR=$PWD
MappedLen=500
MappedIdentity=80
while getopts o:g:a:t:M:I: option
do
        case "${option}"
        in
        a) Fasta=${OPTARG};;
        g) Snova=${OPTARG};;
        o) DIR=${OPTARG};;
	t) Thread=${OPTARG};;
	M) MappedLen=${OPTARG};;
        I) MappedIdentity=${OPTARG};;
        esac
done

if [ -z "$Fasta" ] || [ -z "$Snova" ]; then
        echo "
	usage: /opt/GCB_GapFilling/Fill.sh -a FASTA -g G_CONTIGS [-t THREADS] -o OUTDIR
       	[-M MIN_MAPLEN] [-I MIN_MAPIDENTITY]
 
        positional arguments:
                   -a FASTA
                      Draft assembly FASTA file to be filled
                   -g G_CONTIGS
                      Scaffolds/contigs from other assemblies or long reads for filling in gaps
                   -o OUTDIR
                      Output directory
 
        optional arguments:
                   -t THREADS 
                      Number of threads [default=16]
                   -M MIN_MAPLEN
                      Minimum mapped length of contigs [default=500]
                   -I MIN_MAPIDENTITY
                      Minimum mapping identity of contigs [default=80]
	" 1>&2
	#"Usage: $0 [-a: Fasta to be filled by supernova's contigs] [-s: Supernova fasta file] [-o: Output directory] [-t: Number of threads] [-M: Minimum mapped length of contigs] [-I: Minimum mapping identity of contigs]" 1>&2
        exit 1
fi

SECONDS=0
Rec=$DIR/GCBFilled_Record.log

#Split G-contigs by gaps
echo "[00:00:00]" Removing Ns in input g-contigs > $Rec
python /opt/GCB_GapFilling/RenameSplitbyN.py $Snova
snova=`echo $Snova | awk '{split($1,s,".fa"); print s[1]"_SplitbyN.fa"}'`

#Rename scaffold name from reference
printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Renaming fasta..." >> $Rec
FArename=`echo $(basename "$Fasta") | awk '{split($1,s,".fa"); print s[1]}'`
awk '{if(substr($1,1,1)==">")n=$1;else len[n]+=length($1);}END{for(i in len)print i, len[i];}' $Fasta | sort -k2,2nr | awk '{i++; print $1"\t"">scaffold"i"|size"$2}' > $DIR/${FArename}_RenameBylength.tsv
awk -v dir=$DIR 'NR==FNR{scaf[$1]=$2; next} {if(substr($1,1,1)==">"){print scaf[$1]; split(scaf[$1],t,"|"); sidx=substr(t[1],2); p=scaf[$1];}else{print; p=$1} print p > dir"/"sidx".fa";}' $DIR/${FArename}_RenameBylength.tsv $Fasta > $DIR/${FArename}_rename.fasta

#Mapping
printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Aligning g-contigs to genome..." >> $Rec
bwa index $DIR/${FArename}_rename.fasta
bwa mem -t $Thread -a $DIR/${FArename}_rename.fasta $snova > $DIR/aln.sam

printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Mapping is done. Start to fill the gaps" >> $Rec
out=$DIR"/${FArename}_GCBFilled.fasta"

#split sam to each scaffold's contig.txt
awk -v d=$DIR 'and($2,4)==0 && $1!="@SQ" && $1!="@PG" {split($3,name,"|"); if(and($2,256)==0){  \
	r[$1"_"$3]++; if(r[$1"_"$3]==1) {if($16 && substr($16,1,5)=="SA:Z:"){outSA="SA:Z:"; split(substr($16,6),sa,";"); for(i in sa){split(sa[i],s,","); if(s[1]==$3)outSA=outSA""sa[i]";"}} if(outSA=="SA:Z:")outSA=""; print $1, $2, $4, $5, $6, outSA >> d"/"name[1]"_primary_contig.txt";}} \
	else{print $1, $2, $3, $4, $5, $6, $16 >> d"/"name[1]"_secondary_contig.txt";}}' $DIR/aln.sam
awk '{split($2,s,"|"); print substr(s[1],2)}' $DIR/${FArename}_RenameBylength.tsv > $DIR/ScafName.txt

Filling () {
	DIR=$1
	ScafName=$2
	snova=$3
	MappedLen=$4
       	MappedIdentity=$5

	fa=$DIR/$ScafName.fa
	python /opt/LCB_GapFilling/locateGap.py $fa > $DIR/gap_pos_$ScafName.txt
	gappos=$DIR/gap_pos_$ScafName.txt

	if [ ! -f $DIR/${ScafName}_primary_contig.txt ]; then
		echo $fa >> $DIR/getFA.txt
		rm $gappos
		return
	fi
	if [ -f $DIR/${ScafName}_secondary_contig.txt ]; then
		awk 'NR==FNR{if(and($2,16)==16) f="-"; else f="+"; aln[$1]=$3","$4","f","$6","$5",0;"aln[$1]; next} {if(aln[$1]){if(!$6 || substr($6,1,2)!="SA")$6="SA:Z:"aln[$1]; else $6=$6""aln[$1]} print $0;}' $DIR/${ScafName}_secondary_contig.txt $DIR/${ScafName}_primary_contig.txt >  $DIR/${ScafName}_contig.txt
		rm $DIR/${ScafName}_secondary_contig.txt $DIR/${ScafName}_primary_contig.txt
	else
		mv $DIR/${ScafName}_primary_contig.txt $DIR/${ScafName}_contig.txt
	fi
	python /opt/LCB_GapFilling/contigRange.py $DIR/${ScafName}_contig.txt > $DIR/${ScafName}_contigRange.txt
	gapCnt=`wc -l $gappos | awk '{print $1;}'`
	if [ $gapCnt -eq "0" ]; then
		echo $fa >> $DIR/getFA.txt
		rm $gappos
		return
	fi

	awk 'NR==FNR{gap[i]=$1;gapEnd[i]=$2;gapC[i]="{";gapPL[i]="{";gapPR[i]="{";i++;next} {s[$1]["complete"]="{"; s[$1]["partial"]="{"; for(i in gap){ \
	if(gapEnd[i]<$4 && gap[i]>$3){ s[$1]["complete"]=s[$1]["complete"]""gap[i]";"; gapC[i]=gapC[i]""$1";";} \
	else if((gap[i]-$4<0?$4-gap[i]:gap[i]-$4)<=50){s[$1]["partial"]=s[$1]["partial"]""gap[i]";"; gapPL[i]=gapPL[i]""$1";";} \
	else if ((gapEnd[i]-$3<0?$3-gapEnd[i]:gapEnd[i]-$3)<=50){s[$1]["partial"]=s[$1]["partial"]""gap[i]";"; gapPR[i]=gapPR[i]""$1";";} } \
	}END{for(j in gap) print gap[j]" "gapEnd[j]" complete: "gapC[j]"} partialL: "gapPL[j]"} partialR: "gapPR[j]"}"}' \
	$gappos $DIR/${ScafName}_contigRange.txt | sort -k1,1n > $DIR/${ScafName}_coveredContig.txt
		
	python /opt/LCB_GapFilling/filling.py $DIR/${ScafName}_coveredContig.txt $DIR/${ScafName}_contigRange.txt $fa $snova $MappedLen $MappedIdentity > $DIR/${ScafName}_filling.log 
        mv $DIR/${ScafName}_filled_fix.fa $DIR/${ScafName}_filled.fa
	echo $DIR/${ScafName}_filled.fa >> $DIR/getFA.txt

	python /opt/LCB_GapFilling/locateGap.py  $DIR/${ScafName}_filled.fa > $DIR/gap_pos_${ScafName}_filled.txt
	Total_gaps=$(awk '{print $1;}' $gappos | wc -l)
	echo 'Total gaps:' $Total_gaps > $DIR/${ScafName}_gapFilling.log

	R=$DIR/gapStatus_Record.txt
	Filled_gaps=$(awk '$6=="Filled"{print $1;}' $R | uniq | wc -l)
	F_per=$(echo $Filled_gaps $Total_gaps | awk '{printf("%.2f", ($1/$2)*100)}')
	UF_per=$(echo $F_per | awk '{printf("%.2f", 100-$1)}')

	echo 'Filled gaps:' $Filled_gaps '('$F_per'%)' >> $DIR/${ScafName}_gapFilling.log
	echo '---Completely filled gaps:' $(awk '$3=="Complete"{print $1;}' $R | wc -l) >> $DIR/${ScafName}_gapFilling.log
	echo '---partially filled gaps:' $(awk '$3!="Complete" && $3!="X" && $3!="Delete"{print $1;}' $R | uniq | wc -l) >> $DIR/${ScafName}_gapFilling.log
	echo '---Delete gaps:' $(awk '$3=="Delete"{print $1;}' $R | uniq | wc -l) >> $DIR/${ScafName}_gapFilling.log
	echo ' ' >> $DIR/${ScafName}_gapFilling.log
	echo 'Unfilled gaps:' $((Total_gaps-Filled_gaps)) '('$UF_per'%)' >> $DIR/${ScafName}_gapFilling.log
	echo '---Inconsistent gaps:' $(awk '$6=="Inconsistent" || $6=="inDeletion"{print $1;}' $R | wc -l) >> $DIR/${ScafName}_gapFilling.log
	echo '---Conflict gaps:' $(awk '$6=="Conflict"{print $1;}' $R | wc -l) >> $DIR/${ScafName}_gapFilling.log
	echo '---Unfillale gaps:' $(awk '$6=="noContig"{print $1;}' $R | wc -l) >> $DIR/${ScafName}_gapFilling.log
	echo ' ' >> $DIR/${ScafName}_gapFilling.log
	Total_gapLength=$(awk '{if($4=="")len+=$3; else len+=$4;}END{print len;}' $gappos)
	echo "Total Gap Length in Scaffold: " $Total_gapLength " (bp)" >> $DIR/${ScafName}_gapFilling.log
	Remain_gapLength=$(awk 'BEGIN{len=0} {len+=$3}END{print len;}' $DIR/gap_pos_${ScafName}_filled.txt)
	echo "Remaining Gap Length After Gap-Filling: " $Remain_gapLength " (bp)" >> $DIR/${ScafName}_gapFilling.log
	echo ">>> " $(echo $Total_gapLength $Remain_gapLength | awk '{printf("%.2f", ($1-$2)*100/$1)}') "% of Gaps are filled! <<<" >> $DIR/${ScafName}_gapFilling.log
	echo $ScafName is done!
	rm $fa $gappos $DIR/${ScafName}_coveredContig.txt $DIR/gap_pos_${ScafName}_filled.txt $DIR/${ScafName}_filling.log $DIR/${ScafName}_contig.txt
}
export -f Filling

cat $DIR/ScafName.txt | xargs -n 1 -I {} -P $Thread bash -c 'Filling $1 $2 $3 $4 $5' sh $DIR {} $snova $MappedLen $MappedIdentity
cat $DIR/getFA.txt | xargs -n 1 cat > $out
rm $DIR/*_secondary_contig.txt $DIR/scaffold*.fa $DIR/getFA.txt $DIR/aln.sam $DIR/ScafName.txt

printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Gap filling with g-contigs is done!" >> $Rec
echo Output fasta file is at $out >> $Rec
