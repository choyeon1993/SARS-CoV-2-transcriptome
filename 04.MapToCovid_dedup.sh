directory=/scratch/2020-06-04/bio-zhaoy1/Ribozero/3.star_mapping_to_covid
outdir=/scratch/2020-06-04/bio-zhaoy1/Ribozero/3.star_mapping_to_covid/dedup
mkdir -p $outdir
cd $directory
DIR=/scratch/2020-06-04/bio-zhaoy1/90.scripts/star_jobs
samples=`ls V*_Aligned.sortedByCoord.out.bam`
for sample in $samples
do
samplename=`basename $sample _Aligned.sortedByCoord.out.bam`
#fq1=$sample
#fq2=$samplename\_Unmapped.out.mate2.fastq.gz
cd $DIR
config=$samplename\_ribozeroToVirus_dedup.lsf
touch $config

echo "#extract unmapped reads from map-to-vero bam files and map the reads to virus genome" > $config
echo "#!/bin/bash" > $config
#echo "#BSUB -J "$samplename >> $config
echo "#BSUB -q short" >> $config
echo "#BSUB -n 40" >> $config
echo "#BSUB -R \"span[hosts=1]\"" >> $config
echo "#BSUB -W 120:00" >> $config
echo "#BSUB -o "$samplename"_riboToVirus_dedup.out" >> $config
echo "#BSUB -e "$samplename"_riboToVirus_dedup.err" >> $config

echo "directory="$directory >> $config
echo "sample="$sample >> $config
echo "outdir="$outdir >> $config
echo "samplename="$samplename >> $config
echo "base="$samplename >> $config
#echo "fq1="$fq1 >> $config
#echo "fq2="$fq2 >> $config
echo "cd \$directory" >> $config
echo "echo \$directory/\$sample" >> $config
echo "samtools view -bf 4 \$sample > \$samplename.unmapped.bam" >> $config
echo "outdir="$outdir >> $config
echo "cd \$outdir" >> $config
echo "fqdir=\$outdir/fqfiles" >> $config
echo "mkdir -p \$fqdir" >> $config
echo "bamToFastq -i \$directory/\$samplename.unmapped.bam -fq \$fqdir/\$samplename.unmapped.1.fq -fq2 \$fqdir/\$samplename.unmapped.2.fq" >> $config
echo "gzip \$fqdir/\$samplename.unmapped*" >> $config
echo "genomedir=/work/bio-zhaoy1/03.Reference/COVID19/STAR_index" >> $config

echo "#mark duplicated reads" >> $config
echo "java -jar /work/bio-zhaoy1/90.software/picard.jar MarkDuplicates I=\$sample O=\$outdir/\$base.marked_dup.bam M=\$outdir/\$base.dup_metrics.txt" >> $config
echo "#filter duplicated reads" >> $config
echo "samtools view -b -@ 8 -o \$outdir/\$base.no_dup.bam -F 1024 \$outdir/\$base.marked_dup.bam" >> $config

echo "cd "$outdir >> $config
echo "jumpdir=jumps" >> $config
echo "mkdir -p \$jumpdir" >> $config

echo "#generate index files" >> $config
echo "samtools index \$base.no_dup.bam > \$base.no_dup.bam.bai" >> $config
<<!
#for Poly(A)-seq data
echo "#extract cigar containing jumps" >> $config
echo "samtools view \$base.no_dup.bam | awk 'BEGIN {{ OFS=\"\t\"; }} {{ if (\$6 ~ /N/) print \$4,\$6; }}' > \$jumpdir/\$base.cigar_jumps.txt" >> $config

echo "#counting jump sites from Poly(A) Sequencing data" >> $config
echo "cat \$jumpdir/\$base.cigar_jumps.txt |python /scratch/2020-06-04/bio-zhaoy1/90.scripts/convert-cigars.py| sort -k1,2n | uniq -c | awk '{{ OFS=\"\t\"; }} {{ print \$2,\$3,\$1; }}' > \$jumpdir/\$base.jumps.txt" >> $config
!
echo "#calculate coverage" >> $config
echo "bedtools genomecov -ibam  \$base.no_dup.bam -dz -split > \$base.dedup.mapped.coverage.txt" >> $config
echo "#pileup reads" >> $config
echo "genomefa=/work/bio-zhaoy1/03.Reference/COVID19/COVID19.fasta" >> $config
echo "samtools mpileup -f \$genomefa \$base.no_dup.bam > \$base.dedup.sort.mapped.pileup.txt" >> $config

echo "#extract junction-only reads from bam" >> $config
echo "(samtools view -H \$base.no_dup.bam; samtools view -q 100 \$base.no_dup.bam |awk '{{ if (\$6 ~ /N/) print \$0; }}') > \$base.temp.bam" >> $config
echo "samtools view -f \$base.temp.bam -o \$base.sort.dedup.mapped.junctionOnly.bam" >> $config
echo "rm \$base.temp.bam" >> $config

echo "cd \$outdir" >> $config
echo "bedtools genomecov -ibam \$base.no_dup.bam -dz -split -du -strand + > \$base.dedup.mapped.coveragebyStrand.txt" >> $config
echo "bedtools genomecov -ibam \$base.no_dup.bam -dz -split -du -strand - > \$base.dedup.mapped.coveragebyNegativeStrand.txt" >> $config

echo "cd "$outdir >> $config
echo "jumpdir=jumps-bystrand" >> $config
echo "mkdir -p \$jumpdir" >> $config

#####for Ribozero seq data
echo "#extract cigar containing jumps" >> $config
echo "samtools view \$base.no_dup.bam | awk 'BEGIN {{ OFS=\"\t\"; }} {{ if (\$6 ~ /N/) print \$2,\$4,\$6; }}' > \$jumpdir/\$base.cigar_jumps.txt" >> $config

echo "#counting jump sites from Ribozero data" >> $config
echo "cat \$jumpdir/\$base.cigar_jumps.txt |python /scratch/2020-06-04/bio-zhaoy1/90.scripts/convert-cigars-singlestrand.py| sort -k2,3n | uniq -c | awk '{{ OFS=\"\t\"; }} {{ print \$3,\$4,\$2,\$1; }}' > \$jumpdir/\$base.jumps.txt" >> $config


#echo "jumpdir=no_dup_jumps" >> $config
#echo "mkdir -p $jumpdir" >> $config
#echo "#extract cigar containing jumps from no_duplication bam file" >> $config
#echo "samtools view \$base.no_dup.bam | awk \'BEGIN {{ OFS=\"\t\"; }} {{ if (\$6 ~ /N/) print \$4,\$6; }}\'> \$jumpdir/\$base.cigar_jumps_no_dup.txt" >> $config
#echo "#counting jump sites" >> $config
#echo "cat \$jumpdir/\$base.cigar_jumps_no_dup.txt |python /home/bio-zhaoy1/03.Projects/04.COVID/90.scripts/convert-cigars.py| sort -k1,2n | uniq -c | awk \'{{ OFS=\"\t\"; }} {{ print \$2,\$3,\$1; }}\' > \$jumpdir/\$base.jumps_no_dup.txt" >> $config


bsub -J $samplename < $config

done
