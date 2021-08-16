indir=/scratch/Shares/public/sread2021/data_files/day9/ATACbams/
outdir=/scratch/Users/maallen3/ATAC/

mkdir -p $outdir

#for pathandfilename in `ls ${indir}*sorted.bam`; do
for pathandfilename in `ls ${indir}*D8E_B*sorted.bam`; do
rootname=`basename $pathandfilename .sorted.bam`
echo $rootname

sbatch --export=indir=$indir,rootname=$rootname,outdir=$outdir HMMRATACandsort.sh
done
