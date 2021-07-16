#!/bin/bash

echo "Process read1 to generate index file..."

create_index="$(gunzip -c "$1" | awk '{if($1~/^@/&&$2~/^1:N:0:/){split($2,si,":");qstr=si[4];if(match(si[4],"[N|n]")){gsub("[N|n]","#",qstr);gsub("[A|T|G|C|a|t|g|c]","F",qstr);print $0"\n"si[4]"\n+\n"qstr}else{gsub("[A|T|G|C|a|t|g|c]","F",qstr);print $0"\n"si[4]"\n+\n"qstr}}}' | gzip > "$2")"
#echo "$create_index"
