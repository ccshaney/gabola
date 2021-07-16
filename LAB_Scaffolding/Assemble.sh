DIR=$PWD
Thread=8
n=2000

while getopts f:o:r:t:n: option
do
        case "${option}" in
                f) BXCntFile=${OPTARG};;
                o) DIR=${OPTARG};;
                r) Read=${OPTARG};;
		t) Thread=${OPTARG};;
		n) Num=${OPTARG};;
        esac
done

if [ -z "$BXCntFile" ] || [ -z "$Read" ]; then
	echo "
	usage: /opt/LAB_Scaffolding/Assemble.sh
               -f SCAFFOLD_PAIR -r FASTQ [-t THREADS] [-n NUM] -o OUTDIR

        positional arguments:
               -f SCAFFOLD_PAIR
                  Candidate scaffold pairs and shared barcode counts
                  (C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv)
               -r FASTQ
                  Directory for non-dup split fastq from Preprocess module step I
               -o OUTDIR
                  Output directory (same as Preprocessing.sh)

        optional arguments:
               -t THREADS
                  Number of tasks/scaffold pairs per run (one task occupies 10 threads);
		  depending on how many available CPUs [default=8]
               -n NUM
                  Process top n scaffold end pairs [default=2000]
       " 1>&2
	#echo "Usage: $0 [-f: ScaffoldEndPair.tsv] [-r: Directory for splited barcoded fastq] [-o: Output directory for assembled contigs] [-t: Number of Jobs per run(each job runs with 10 threads)] [-n: Maximum number of scaffold end pairs to be processed]" 1>&2 
	exit 1
fi

CollectBXReads () {
        if [ ! -f $2/BXList_$1 ]; then
                return
        fi
        while read -r BX; do
                f1=$3/$BX
		cat ${f1}_1 >> $2/R1_$1.fq
                cat ${f1}_2 >> $2/R2_$1.fq

        done < $2/BXList_$1
}
export -f CollectBXReads

assemble () {
	#Directory=$1
	#SpltedReadDirectory=$2
	#scaf1=$3
	#scaf2=$4
	#SharedBXCnt=$5

	s1=`echo $3 | awk '{split($1,a,"_"); print a[1];}'`
        p1=`echo $3 | awk '{print substr($1,length($1)-3,4)}'`
        s2=`echo $4 | awk '{split($1,a,"_"); print a[1];}'`
        p2=`echo $4 | awk '{print substr($1,length($1)-3,4)}'`

	wd=$1/${s1}${p1}_$s2${p2}

	if [ ! -d $wd/ ] && [ ! -d $1/$s2${p2}_$s1$p1/ ]; then
		SECONDS=0
		mkdir $wd
		cat $1/HeadTail_BXList/${3}_BXList $1/HeadTail_BXList/${4}_BXList | sort | uniq > $wd/BXList

		collecting_thread=6
        	awk -v Thread=$collecting_thread -v p=$wd '{print > p"/BXList_"(NR%Thread+1)}' $wd/BXList
	        seq 1 $collecting_thread | xargs -n 1 -I {} -P ${collecting_thread} bash -c 'CollectBXReads $1 $2 $3' sh {} $wd $2
		ls $wd/R1_*.fq | xargs -n 1 cat > $wd/R1.fq
		ls $wd/R2_*.fq | xargs -n 1 cat > $wd/R2.fq
		rm $wd/BXList_*

		BXCnt=$(wc -l $wd/BXList | awk '{print $1;}')
                RdCnt=$(wc -l $wd/R1.fq | awk '{print $1/4;}')

                python /usr/local/bin/SPAdes-3.15.2-Linux/bin/spades.py -t 10 -1 $wd/R1.fq -2 $wd/R2.fq -o $wd/ &> $wd/a.log
                rm $wd/a.log $wd/*.fq $wd/assembly_* $wd/contigs.* $wd/before_rr.fasta
		rm -r $wd/corrected $wd/K*/ $wd/misc

                diff=$SECONDS

                echo ${s1}${p1}_$s2$p2 Done! "#BX:" $BXCnt "#Rp:" $RdCnt Elapsed Time: $(($diff / 3600)):$((($diff / 60) % 60)):$(($diff % 60)) >> $1/Assemble.log
	fi
}

export -f assemble

echo Assembling of $BXCntFile starts with $Thread jobs per run: `date '+%Y-%m-%d %H:%M:%S'` > $DIR/Assemble.log
head -n $n $BXCntFile | xargs -n 1 -I {} -P $Thread bash -c 'assemble $1 $2 $3' sh $DIR $Read {}
echo Finished: `date '+%Y-%m-%d %H:%M:%S'` >> $DIR/Assemble.log
echo Next step: Scaffolding.sh with -f $BXCntFile -o $DIR and -a FASTA >> $DIR/Assemble.log
