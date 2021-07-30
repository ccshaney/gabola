DIR=$PWD
Top=2
while getopts f:o:p:v: option
do
	case "${option}"
	in
	f) BXCntFile=${OPTARG};;
	p) PairSum=${OPTARG};;
	o) DIR=${OPTARG};;
	v) Top=${OPTARG};;
	esac
done

if [ -z "$BXCntFile" ] || [ -z "$PairSum" ]; then
        echo "
	usage: /opt/GCB_Scaffolding/CandidatePair.sh
               -f SCAFFOLD_PAIR -p BXLIST [-v PAIR_NUM] -o OUTDIR

        positional arguments:
               -f SCAFFOLD_PAIR
                  Scaffold pairs and shared barcode counts from Preprocess module step III
		  (C70M60_ScafA_ScafB_BXCnt.tsv)
               -p BXLIST
                  Barcodes on every scaffold end from Preprocess module step III
		  (C70M60_ScafHeadTail_BX_pairSum.tsv)
               -o OUTDIR
                  Output directory

        optional arguments:
               -v PAIR_NUM
                  Remove multiple ends, keep top v pairs [default=2]
        " 1>&2
	#echo "Usage: $0 [-f: XXX_C70M60_ScafA_ScafB_BXCnt.tsv] [-p: XXX_C70M60_ScafHeadTail_BX_pairSum.tsv] [-v: Keep top v scaffold ends] [-o: Output directory for assembled contigs]" 1>&2
        exit 1
fi

#Delete pairs with same scaffold (ex: s1_Head s1_Tail) and  Remove multiple ends, keep top 2 pairs
awk -v t=$Top '{split($1,s1,"_"); split($2,s2,"_"); if(s1[1]!=s2[1]){ if(s[$1]<t && s[$2]<t){print; s[$1]++; s[$2]++;}}}' $BXCntFile > $DIR/C70M60_ScafA_ScafB_BXCnt_rmMultiEnd.tsv

#Create BX list for each Head/Tail
mkdir $DIR/HeadTail_BXList
awk -v dir=$DIR '{print $1 > dir"/HeadTail_BXList/"$2"_BXList";}' $PairSum


