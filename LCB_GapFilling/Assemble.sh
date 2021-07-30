#!/bin/bash
DIR=$PWD
Thread=8
JQ=1

while getopts r:o:t:q: option
do
        case "${option}"
        in
        r) Read=${OPTARG};;
        o) DIR=${OPTARG};;
        t) Thread=${OPTARG};;
	q) JQ=${OPTARG};;
        esac
done

if [ -z "$Read" ]; then
	echo "
	Usage: /opt/LCB_GapFilling/Assemble.sh -r FASTQ [-t THREADS] [-q JobQueue_NUM ] -o OUTDIR
 
              positional arguments:
                         -r FASTQ
                            Directory for non-duplicated split fastq from Preprocess step I
                         -o OUTDIR	
                            Output directory (same as outdir of ProduceBXList.sh)
 
              optional arguments:
                         -t THREADS
                            Number of tasks/gaps per run (one task occupies 10 threads);
			    depending on how many available CPUs [default=8]
                         -q JobQueue_NUM
                            Index of Task queue [default=1]
" 1>&2
	#echo "Usage: $0 [-r: Directory for splited nondup fastq] [-q: Index of job queue] [-o: Output directory for assembled contigs] [-t: Number of Jobs per run(each job runs with 10 threads)]" 1>&2 
	exit 1
fi


CollectBXReads () {
        if [ ! -f $2/BXList_$1 ]; then
                return
        fi
        while read -r BX; do
                #For eel
		#f1=$(echo $BX $3 | awk '{print $2"/"substr($1,1,4)"/"substr($1,5,4)"/"$1}')

		f1=$3/$BX
		cat ${f1}_1 >> $2/R1_$1.fq
                cat ${f1}_2 >> $2/R2_$1.fq

        done < $2/BXList_$1
}
export -f CollectBXReads

assemble () {
	#Directory=$1
	#Splited nondup fastq=$2
	#scafID=$3
	#gapID=$4

	wd=$1/$3/$4
	if [ -d $wd ]; then
		collecting_thread=8
        	awk -v Thread=$collecting_thread -v p=$wd '{print > p"/BXList_"(NR%Thread+1)}' $wd/BXList
	        seq 1 $collecting_thread | xargs -n 1 -I {} -P ${collecting_thread} bash -c 'CollectBXReads $1 $2 $3' sh {} $wd $2
		ls $wd/R1_*.fq | xargs -n 1 cat > $wd/R1.fq
		ls $wd/R2_*.fq | xargs -n 1 cat > $wd/R2.fq
		rm $wd/BXList_* $wd/R1_*.fq $wd/R2_*.fq

		BXCnt=$(wc -l $wd/BXList | awk '{print $1;}')
                RdCnt=$(wc -l $wd/R1.fq | awk '{print $1/4;}')

                python /usr/local/bin/SPAdes-3.15.2-Linux/bin/spades.py -t 10 -1 $wd/R1.fq -2 $wd/R2.fq -o $wd/ &> $wd/a.log
		if [ -f $wd/scaffolds.fasta ]; then
			mv $wd/scaffolds.fasta $1/$3/scaffolds_$4.fasta
			rm -r $wd
		fi

                echo $3/$4 Done! "#BX:" $BXCnt "#Rp:" $RdCnt >> $1/Assemble_Job${5}_Record.log
	fi
}
export -f assemble

echo Assembling starts with $Thread jobs per run: `date '+%Y-%m-%d %H:%M:%S'` > $DIR/Assemble_Job${JQ}_Record.log

if [ ! -f $DIR/JobList$JQ ]; then
	echo "Error: $DIR/JobList$JQ: No such file for following assembling jobs."
	exit 1
else
	while read -r scaf; do
		t=`wc -l $DIR/$scaf/gapRange_BXs.txt | awk '{print $1;}'`
		for ((i=1;i<=t;i++));
		do
			echo $scaf gapID$i >> $DIR/JobList${JQ}_Works
		done	
	done < $DIR/JobList$JQ
fi

cat $DIR/JobList${JQ}_Works | xargs -n 1 -I {} -P $Thread bash -c 'assemble $1 $2 $3 $4 $5' sh $DIR $Read {} $JQ

echo All done: `date '+%Y-%m-%d %H:%M:%S'` >> $DIR/Assemble_Job${JQ}_Record.log
