
#!/bin/bash

#went to https://trace.ncbi.nlm.nih.gov/Traces/study/?go=home and searched for SRP002796 to download both the SRR_ACC_List.txt and the Sra-run-table
#to run type bash downloadall.sh <path_to_SRR_ACC_List.txt> <outdir> <email>

mkdir -p $2

IFS=''
while read var
do
mkdir -p ${2}/eando/
echo $var
if [ -n "$var" ];
then
sbatch -J $var --mail-user=$3 --output=${2}/eando/${var}.%j.out --error=${2}/eando/${var}.%j.err --export=filename=$var,outdir=$2 downloadafastq.sh
fi
done < $1



