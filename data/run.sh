#!/bin/bash

if [ $# -eq 0 ] ; then
    echo "Usage: run.sh <MODE>"
    echo "MODE = conda | docker | wdl"
    exit 1
fi
MODE=$1

#
# run varpipe
#
#FASTQ=$HOME/Analysis/wgs_pipeline/data/
FASTQ=$PWD
R1=$FASTQ/ERR552797_1.fastq.gz
R2=$FASTQ/ERR552797_2.fastq.gz
R1=$FASTQ/ERR552797_30percent_1.fq.gz
R2=$FASTQ/ERR552797_30percent_2.fq.gz
SAMPL=ERR552797
REF=../tools/NC_000962.3.fa
NAME="NC_000962.3"

if [ "$MODE" = "conda" ] ; then
   ../tools/Varpipeline -q $R1 -q2 $R2 -r $REF -g $NAME -n $SAMPL -o $SAMPL -a -v -t 8 -k
fi

#
# run Docker varpipe
#

#TOOLS=/NCHHSTP-DTBE-Varpipe-WGS/tools
TOOLS=$HOME/Analysis/varpipe3/tools
TLS=/mnt/tools
REF=$TLS/NC_000962.3.fa

FASTQ=$HOME/Analysis/wgs_pipeline/data/
FSQ=/mnt/fastq
R1=$FSQ/ERR552797_1.fastq.gz
R2=$FSQ/ERR552797_2.fastq.gz
R1=$FSQ/ERR552797_30percent_1.fq.gz
R2=$FSQ/ERR552797_30percent_2.fq.gz

DATA=$HOME/Analysis/varpipe3/data
WRK=/mnt/data

if [ "$MODE" = "docker" ] ; then
   docker run --rm -v $FASTQ:$FSQ -v $TOOLS:$TLS -v $DATA:$WRK -w $WRK dbest/varpipe3:latest Varpipeline -q $R1 -q2 $R2 -r $REF -g $NAME -n $SAMPL -v -a -t 8 -k
fi

#
# run WDL script varpipe
#
if [ "$MODE" = "wdl" ] ; then
   miniwdl run varpipe.wdl -i varpipe.yaml
fi

exit 0

