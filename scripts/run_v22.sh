set -euo pipefail
W=/media/nighthawk/_mnt_16TB/webtool/validation
IDX=/media/nighthawk/_mnt_3TB/ref_genomes/GRCh38_STAR-2.7.10b_GTF101_oh100
V22=$W/ref/gencode.v22.annotation.gtf
cd $W
mkdir -p bam22 counts logs
for S in HepG2 sigmoid_colon; do
  B=bam22/${S}_v22_Aligned.out.bam
  if [ ! -s "$B" ]; then
    echo "[$(date +%H:%M)] STAR(v22) $S"
    STAR --genomeDir "$IDX" --sjdbGTFfile "$V22" --sjdbOverhang 100 \
         --readFilesIn fastq/${S}_r1.fastq.gz fastq/${S}_r2.fastq.gz --readFilesCommand zcat \
         --runThreadN 48 --outSAMtype BAM Unsorted --outSAMattributes NH HI AS NM MD \
         --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 \
         --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
         --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --quantMode GeneCounts \
         --limitSjdbInsertNsj 2000000 --outFileNamePrefix bam22/${S}_v22_ > logs/star22_${S}.log 2>&1
  fi
  echo "[$(date +%H:%M)] htseq ine v22-on-v22 $S"
  [ -s counts/${S}_htseq_ine_v22on22.tsv ] || htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id \
      -m intersection-nonempty -n 16 "$B" "$V22" > counts/${S}_htseq_ine_v22on22.tsv 2> logs/htseq_v22on22_${S}.log
done
echo "[$(date +%H:%M)] v22 arm done"
