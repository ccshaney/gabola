#!/bin/bash

echo "Process read1 to generate index file..."

create_index_noZip="$(awk '{if($1~/^@/){split($2,si,":");qstr=si[4];if(match(si[4],"[N|n]")){gsub("[N|n]","#",qstr);gsub("[A|T|G|C|a|t|g|c]","F",qstr);print $0"\n"si[4]"\n+\n"qstr}else{gsub("[A|T|G|C|a|t|g|c]","F",qstr);print $0"\n"si[4]"\n+\n"qstr}}}' $1 | gzip > $2)"
#echo "$create_index_noZip"
