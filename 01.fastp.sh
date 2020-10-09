
directories=(/scratch/2020-04-22/bio-zhaoy1/01.Rawdata/X101SC20030251-Z01-F004/00.Rawdata /scratch/2020-04-22/bio-zhaoy1/01.Rawdata/X101SC20030251-Z01-F005/00.Rawdata /scratch/2020-04-22/bio-zhaoy1/01.Rawdata/X101SC20030251-Z01-J002/00.Raw_data)
for directory in ${directories[@]}
do
#directory=/home/shareDir/bio-share-chenwei/data/Yan/01.COVID19/01.Rawdata/X101SC20030251-Z01-J002/00.Raw_data
#directory=/home/shareDir/bio-share-chenwei/data/Yan/01.COVID19/01.Rawdata/X101SC20030251-Z01-F006/00.Rawdata
cd $directory
samples=`ls -d H*`
for i in $samples
do
#if [ "$i" != "HB2" ]; then
#DIR=/home/shareDir/bio-share-chenwei/data/Yan/01.COVID19/90.scripts
DIR=/scratch/2020-04-22/bio-zhaoy1/90.script/1.fastp_jobs
mkdir -p $DIR
sampledir=$directory/$i
samplename=`basename $sampledir`
echo $samplename
cd $DIR
config=$DIR/$samplename.pbs
touch $config

echo '##fastp deal with raw data' > $config
echo "#PBS -N "$samplename  >> $config
echo "#PBS -l nodes=1:ppn=8"  >> $config
echo "#PBS -j oe"  >> $config
echo "#PBS -o "$samplename"_pbs.log"  >> $config
echo "#PBS -q ser"  >> $config
#echo "#PBS -l nodes=1:ppn=8" >> $config
echo "export  OMP_NUM_THREADS=8" >> $config
echo "#PBS -V"  >> $config
echo "folder="$sampledir >> $config
echo "outdir=/scratch/2020-04-22/bio-zhaoy1/02.cleandata" >> $config
echo "fq1=\$folder/*1.fq.gz" >> $config
echo "fq2=\$folder/*2.fq.gz" >> $config
#echo "base=\`basename \$fq1 _1.fq.gz\`" >> $config
#echo "fastp -n 10 -q 35 -w 8 --detect_adapter_for_pe -i \$fq1 -I \$fq2 -o \$outdir/"$samplename"\_1.clean.fq.gz -O \$outdir/"$samplename"\_2.clean.fq.gz -j \$outdir/"$samplename".json -h \$outdir/"$samplename".html -R "$samplename  >> $config
echo "fastp -n 10 -q 35 -w 8 --overrepresentation_analysis --adapter_fasta /scratch/2020-04-22/bio-zhaoy1/90.script/adapters.fasta -i \$fq1 -I \$fq2 -o \$outdir/"$samplename"\_1.clean.fq.gz -O \$outdir/"$samplename"\_2.clean.fq.gz -j \$outdir/"$samplename".json -h \$outdir/"$samplename".html -R "$samplename  >> $config

qsub $config
#fi
done
#done
