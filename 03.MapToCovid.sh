directory=/scratch/2020-04-26/bio-zhaoy1/2.star_mapping_to_vero
outdir=/scratch/2020-04-26/bio-zhaoy1/3.star_mapping_to_covid
mkdir -p $outdir
cd $directory
DIR=/scratch/2020-04-26/bio-zhaoy1/90.scripts/star_jobs
samples=`ls HA4*.sortedByCoord.out.bam`
for sample in $samples
do
samplename=`basename $sample _Aligned.sortedByCoord.out.bam`

cd $DIR
config=$samplename.lsf
touch $config

#echo "#extract unmapped reads from map-to-vero bam files and map the reads to virus genome" > $config
echo "#!/bin/bash" > $config
#echo "#BSUB -J "$samplename >> $config
echo "#BSUB -q short" >> $config
echo "#BSUB -n 40" >> $config
echo "#BSUB -R \"span[hosts=1]\"" >> $config
echo "#BSUB -W 120:00" >> $config
echo "#BSUB -o "$samplename".out" >> $config
echo "#BSUB -e "$samplename".err" >> $config

echo "directory="$directory >> $config
echo "sample="$sample >> $config
echo "samplename="$samplename >> $config
echo "cd \$directory" >> $config
echo "samtools view -bf 4 \$directory/\$sample > \$directory/\$samplename.unmapped.bam" >> $config
#echo "outdir="$outdir >> $config
#echo "cd \$outdir" >> $config
#echo "fqdir=\$outdir/fqfiles" >> $config
#echo "mkdir -p \$fqdir" >> $config
#echo "bamToFastq -i \$directory/\$samplename.unmapped.bam -fq \$fqdir/\$samplename.unmapped.1.fq -fq2 \$fqdir/\$samplename.unmapped.2.fq" >> $config
#echo "gzip \$fqdir/\$samplename.unmapped*" >> $config
#echo "genomedir=/work/bio-zhaoy1/03.Reference/COVID19/STAR_index" >> $config

#echo "STAR --runThreadN 40 --genomeDir \$genomedir --readFilesCommand zcat --readFilesPrefix \$fqdir/ --readFilesIn \$samplename.unmapped.1.fq.gz \$samplename.unmapped.2.fq.gz --outSAMtype BAM SortedByCoordinate --alignEndsType EndToEnd --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --outSJfilterOverhangMin 12 12 12 12 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterCountTotalMin 1 1 1 1 --outSJfilterDistToOtherSJmin 0 0 0 0 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --scoreGapNoncan -4 --scoreGapATAC -4 --chimOutType WithinBAM HardClip --chimScoreJunctionNonGTAG 0 --alignSJstitchMismatchNmax -1 -1 -1 -1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix "$outdir"/"$samplename"_ToCovid_">> $config

bsub -J $samplename < $config

done
