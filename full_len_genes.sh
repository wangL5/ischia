#!/bin/bash
#pull out full length genes from OG fasta files
#add something to the end of each line
#replace /n with nothing to get everything in one line
#pull out partial=00
#make new lines by replacing the first thing you added to each line with /n.

sed 's/$/:/g' m_S1.fna | awk -v RS="" '{gsub (/\n/,"")}1' | awk -v RS="" '{gsub (/>N/,"\n>N")}1' | grep "partial=00" | sed 's/:/\n/g' | sed '/^$/d' > m_S1.fullgenes.fna
