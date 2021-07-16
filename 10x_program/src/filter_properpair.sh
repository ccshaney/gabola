#!/bin/bash
proper="sambamba view -F 'not (unmapped or mate_is_unmapped or secondary_alignment or failed_quality_control or duplicate or supplementary or chimeric)' -S $1 -o $2"
echo "Start to filter sam file..."
echo "COMMAND: $proper"
comm="$(sambamba view -F "not (unmapped or mate_is_unmapped or secondary_alignment or failed_quality_control or duplicate or supplementary or chimeric)" -S $1 -o $2)"
echo "$comm"
