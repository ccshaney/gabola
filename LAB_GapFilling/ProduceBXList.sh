#!/bin/bash
DIR=$PWD
Rp=3
Group=8
Thread=16
ListSize=10
FlankingRp=2

while getopts a:f:o:n:rp:b:q:t:c:s: option
do
	case "${option}"
	in
	a) Fasta=${OPTARG};;
	f) SamFile=${OPTARG};;
	o) DIR=${OPTARG};;
	n) Num=${OPTARG};;
	rp) Rp=${OPTARG};;
	b) BarcodeNonDupReadPairCnt=${OPTARG};;
	q) Group=${OPTARG};;
	t) Thread=${OPTARG};;
	c) ListSize=${OPTARG};;
	s) FlankingRp=${OPTARG};;
	esac
done

if [ -z "$Fasta" ] || [ -z "$SamFile" ]; then
        echo "
	Usage: /opt/LAB_GapFilling/ProduceBXList.sh -f SAMFILE -a FASTA -o OUTDIR [-n NUM ] 
	[-t THREADS] [-rp MIN_READPAIR_onScaf] [-b BarcodeNonDupReadPairCnt.txt ] [-q  JobQueue] 
	[-c MIN_BarcodeList ] [-s MIN_READPAIR_onGap]
	
	positional arguments:
                -f SAMFILE
                   High quality read-to-scaffold alignment produced from Preprocess Step II
		   (C70M60.sam)
                -a FASTA 
                   Draft assembly with gaps to be filled
                -o OUTDIR
                   Output directory
 
       optional arguments:
                -n NUM
                   Only perform Gap Filling on top n scaffolds, [default=ALL]
                -t THREADS
                   Number of threads [default=16]
                -rp MIN_READPAIR_onScaf
                   Minimum number of read pairs for barcode on scaffold [default=3]
                -s MIN_READPAIR_onGap
                   Minimum number of read pair for barcode on gap’s flanking [default=2]
                -c MIN_BarcodeList
                   Minimum size of gap’s barcode list [default=10]
                -b BarcodeNonDupReadPairCnt.txt
                   BarcodeNonDupReadPairCnt.txt by preprocessing of reads;
		   only required if there are more than one job queue
                -q  JobQueue
                    Number of job queues for gap-filling [default=1]; required if you want to 
		    run GABOLA in parallel on several machines.
		    For example, if you set the parameter as “-q 3”,
		    our program will produce three lists: JobList1, JobList2 and JobList3.
		    Each list containing the gaps to be filled with the format of <scafname><gapID>
" 1>&2
	#echo "Usage: $0 [-f: Read-to-scaffold alignmnt] [-a: Fasta for gap-filling] [-n: Precossing only top n scaffolds] [-t: Number of threads] [-rp: Minimum number of read pairs for barcode on scaffold] [-b: BarcodeNonDupReadPairCnt.txt] [-q: Number of job queues for assembling] [-c: Minimum size of gap's barcode list] [-s: Minimum number of read pair for barcode on gap's flanking] [-o: Output directory for assembled contigs]" 1>&2
        exit 1
fi

SECONDS=0
Rec=$DIR/ProduceBXList_Record.log

#Rename scaffolds
FArename=`echo $(basename "$Fasta") | awk '{split($1,s,".fa"); print s[1]}'`
echo "[00:00:00] Renaming $Fasta; Renamed Fasta is @$DIR/${FArename}_rename.fa" > $Rec
awk '{if(substr($1,1,1)==">")n=$1;else len[n]+=length($1);}END{for(i in len)print i, len[i];}' $Fasta | sort -k2,2nr | awk '{i++; print $1"\t"">scaffold"i"|size"$2}' > $DIR/${FArename}_RenameBylength.tsv
awk 'NR==FNR{scaf[$1]=$2; next} {if(substr($1,1,1)==">")print scaf[$1];else print}' $DIR/${FArename}_RenameBylength.tsv $Fasta > $DIR/${FArename}_rename.fa


#Processing all the scaffolds by default
if [ -z "$Num" ]; then
	Num=`wc -l $DIR/${FArename}_RenameBylength.tsv | awk '{print $1}'`
fi


#Extract n scaffolds from fasta
printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Spliting fasta... $Num scaffolds will be processed." >> $Rec
awk -v n=$Num -v dir=$DIR '{if(substr($1,1,1)==">"){split($1,scaf,"|"); idx=int(substr(scaf[1],10));}  if(idx<=n)print $1 > dir"/scaffold"idx".fa";else print $1 > dir"/Reference_rename_afterScaf"n".fa";}' $DIR/${FArename}_rename.fa


#Get C70M60 read for each scaffold
printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Processing $SamFile..." >> $Rec
python /opt/LAB_GapFilling/splitMappedRead_byScaffold.py $DIR/${FArename}_RenameBylength.tsv $SamFile $Num

if [ ! -z "$BarcodeNonDupReadPairCnt" ]; then
	sort -k1,1 $BarcodeNonDupReadPairCnt > $DIR/BarcodeNonDupReadPairCnt_sort
else
	Group=1
fi


ProduceBXList () {
	DIR=$1
	i=$2
	Rp=$3
	minBXListSize=$4
	minFlankingRp=$5
	BarcodeNonDupReadPairCnt=$6

	mkdir $DIR/scaf$i
	mv $DIR/scaffold$i.fa $DIR/scaf$i/
	mv $DIR/scaffold${i}_C70M60ReadPos.tsv $DIR/scaf$i/
	
	python /opt/LAB_GapFilling/locateGap.py $DIR/scaf$i/scaffold$i.fa > $DIR/scaf$i/gap_pos_scaf$i.txt

	awk -v rp=$Rp 'NR==FNR{bx[$1]++; next} bx[$1]>=2*rp{print;}' $DIR/scaf$i/scaffold${i}_C70M60ReadPos.tsv $DIR/scaf$i/scaffold${i}_C70M60ReadPos.tsv | sort -k3,3nr > $DIR/scaf$i/scaffold${i}_rp${Rp}BX_C70M60ReadPos.tsv
	python /opt/LAB_GapFilling/collectBX_byGap.py $DIR/scaf$i $Rp $minBXListSize $minFlankingRp
	if [ -f $DIR/scaf$i/gapRange_BXs.txt ]; then
                while read -r idx g1 g2 r1 r3 BX; do
                        if [ ! -d $DIR/scaf$i/gapID$idx ]; then
        			mkdir $DIR/scaf$i/gapID$idx
                                echo $DIR"/scaf"$i"/gapID"$idx"/BXList" $BX | awk '{split($2, b, ","); for(j in b)print b[j] > $1;}'
                        fi
                done < $DIR/scaf$i/gapRange_BXs.txt

		if [ ! -z "$BarcodeNonDupReadPairCnt" ]; then
	                awk '{split($6,bx,","); for(i in bx)rec[bx[i]]++;}END{for(i in rec)print i"\t"rec[i];}' $DIR/scaf$i/gapRange_BXs.txt | sort -k1,1 > $DIR/scaf$i/TotalProcessedBXCnt.tsv
	        	python /opt/LAB_GapFilling/computeProcessedReadPairCnt.py $DIR/scaf$i/TotalProcessedBXCnt.tsv $DIR/BarcodeNonDupReadPairCnt_sort >> $DIR/Scaf_ProcessedReadPairCnt.tsv
			rm $DIR/scaf$i/TotalProcessedBXCnt.tsv
		fi
        fi
        rm $DIR/scaf$i/scaffold${i}_rp${Rp}BX_C70M60ReadPos.tsv 
}
export -f ProduceBXList


#Get barcode list on every scaffolds' gaps
printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Extracting BX list for each gap on total $Num scaffolds with $Thread threads..." >> $Rec
seq 1 $Num | xargs -n 1 -I {} -P $Thread bash -c 'ProduceBXList $1 $2 $3 $4 $5 $6' sh $DIR {} $Rp $ListSize $FlankingRp $BarcodeNonDupReadPairCnt


printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "Delivering $Num scaffolds to $Group job queue..." >> $Rec
if [ ! -z "$BarcodeNonDupReadPairCnt" ];then
	sort -k2,2nr $DIR/Scaf_ProcessedReadPairCnt.tsv | awk '$2>0' > $DIR/Scaf_ProcessedReadPairCnt_sorted.tsv
	mv $DIR/Scaf_ProcessedReadPairCnt_sorted.tsv $DIR/Scaf_ProcessedReadPairCnt.tsv
	awk -v group=$Group -v dir=$DIR 'BEGIN{for(i=1;i<=group;i++){record[i]=0; scaf[i]="";}} {min=record[1]; g=1; for(i=2;i<=group;i++){if(record[i]<min){min=record[i]; g=i;}} record[g]+=$2; scaf[g]=scaf[g]""$1"\t"; print $1 > dir"/JobList"g;}END{for(i=1;i<=group;i++)print scaf[i]}' $DIR/Scaf_ProcessedReadPairCnt.tsv > $DIR/JobQueueByScaffold.tsv
	rm $DIR/BarcodeNonDupReadPairCnt_sort
else
	echo $Num $DIR | awk '{for(i=1;i<=$1;i++){s=s"""scaf"i"\t"; print "scaf"i > $2"/JobList1";} print s;}' > $DIR/JobQueueByScaffold.tsv
fi

printf '[%02d:%02d:%02d] %s\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)) "All Done! Output directory is @$DIR" >> $Rec
