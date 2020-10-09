
#directories=(/scratch/2020-04-22/bio-zhaoy1/01.Rawdata/X101SC20030251-Z01-F004/00.Rawdata /scratch/2020-04-22/bio-zhaoy1/01.Rawdata/X101SC20030251-Z01-F005/00.Rawdata /scratch/2020-04-22/bio-zhaoy1/01.Rawdata/X101SC20030251-Z01-J002/00.Raw_data)
#for directory in ${directories[@]}
#do
#directory=/home/shareDir/bio-share-chenwei/data/Yan/01.COVID19/01.Rawdata/X101SC20030251-Z01-J002/00.Raw_data
#directory=/home/shareDir/bio-share-chenwei/data/Yan/01.COVID19/01.Rawdata/X101SC20030251-Z01-F006/00.Rawdata
#directory=/scratch/2020-05-26/bio-zhaoy1/04.ribozero_cleandata
directory=/scratch/2020-05-26/bio-zhaoy1/01.cleandata/
outdir=/scratch/2020-05-26/bio-zhaoy1/Ribozero/2.star_mapTovero
mkdir -p $outdir
cd $directory
samples=`ls V*_1.clean.fq.gz`
for i in $samples
do
#if [ "$i" != "HB2" ]; then
#DIR=/home/shareDir/bio-share-chenwei/data/Yan/01.COVID19/90.scripts
DIR=/scratch/2020-05-26/bio-zhaoy1/90.scripts/star_jobs/ribozero
mkdir -p $DIR
samplename=`basename $i _1.clean.fq.gz`
echo $samplename
fq1=$i
fq2=$samplename\_2.clean.fq.gz
cd $DIR

genomedir=/work/bio-zhaoy1/03.Reference/Chlsabl/star_index
bsub -J $samplename -q short -n 40 -o $samplename\_pbs.log STAR --runThreadN 40 --genomeDir $genomedir --readFilesCommand zcat --readFilesIn $directory/$fq1 $directory/$fq2 --outSAMtype BAM SortedByCoordinate --alignEndsType EndToEnd --outReadsUnmapped Fastx --outFileNamePrefix $outdir/$samplename\_

done
