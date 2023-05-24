#!/bin/bash

#
# run varpipe
#
FASTQ=$HOME/Analysis/wgs_pipeline/data/
R1=$FASTQ/ERR552797_1.fastq.gz
R2=$FASTQ/ERR552797_2.fastq.gz
R1=$FASTQ/ERR552797_30percent_1.fq.gz
R2=$FASTQ/ERR552797_30percent_2.fq.gz
SAMPL=ERR552797
REF=../tools/ref2.fa

#../tools/Varpipeline -q $R1 -q2 $R2 -r $REF -g NC_000962 -n $SAMPL -o $SAMPL -a -v -t 8 -k

#
# run Docker varpipe
#

#TOOLS=/NCHHSTP-DTBE-Varpipe-WGS/tools
TOOLS=$HOME/Analysis/varpipe/tools
TLS=/mnt/tools
REF=$TLS/ref2.fa

FASTQ=$HOME/Analysis/wgs_pipeline/data/
FSQ=/mnt/fastq
R1=$FSQ/ERR552797_1.fastq.gz
R2=$FSQ/ERR552797_2.fastq.gz
R1=$FSQ/ERR552797_30percent_1.fq.gz
R2=$FSQ/ERR552797_30percent_2.fq.gz

DATA=$HOME/Analysis/varpipe3/data
WRK=/mnt/data

#docker run --rm -v $FASTQ:$FSQ -v $TOOLS:$TLS -v $DATA:$WRK -w $WRK dbest/varpipe3:latest Varpipeline -q $R1 -q2 $R2 -r $REF -g NC__000962 -n $SAMPL -v -a -t 8 -k

#
# run WDL script varpipe
#
miniwdl run varpipe.wdl -i varpipe.yaml

exit 0

