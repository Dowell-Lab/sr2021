
Vector[3]=file3.txt
Vector[0]=file1.txt
Vector[1]=file2.txt
Vector[2]=file3.txt

for index in $(seq 0 3)
do
rootname=${Vector[$index]}
fastqc ${indir}${rootname}.fastq

echo $index
echo $rootname
done
