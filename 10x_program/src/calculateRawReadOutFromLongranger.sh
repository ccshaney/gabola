#!/bin/bash

echo "Calculate number of reads and barcode from raw FASTQ file which is output from longranger basic..."
echo "Raw read after longranger process - Number of read and barcode:" > $2
Cal_readNumber="$(gunzip -c "$1" | awk 'BEGIN{i=0;}{if($1~/^@/){i+=1;if($2~/BX:Z:/){gsub("BX:Z:","",$2);gsub("-1","",$2);print $2 >> "temp_10X_barcodeOfRawFQ.txt"}}}END{print "Number of read: "i}' >> "$2")"
Cal_barcodeNumber="$(uniq temp_10X_barcodeOfRawFQ.txt | wc -l)"
echo "Number of barcode: $Cal_barcodeNumber\n" >> $2
echo "DONE!"
