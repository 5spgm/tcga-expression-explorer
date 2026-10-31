#!/usr/bin/env bash
# =====================================================================
# Direct test of the counting-rule mechanism on public RNA-seq
# =====================================================================
# One alignment per sample, then the same BAM counted six ways. The only
# thing that changes between conditions is the counting rule (and, for one
# condition, the annotation version) — so any difference in the TREX1 value
# is attributable to the counting step alone.
#
#   1. STAR --quantMode GeneCounts, GENCODE v36      (reproduces GDC v2)
#   2. htseq-count -m intersection-nonempty, v36     (reproduces the GDC v1 rule)
#   3. htseq-count -m intersection-nonempty, v22     (GDC v1 annotation)
#   4. htseq-count -m union, v36                     (the other HTSeq mode)
#   5. featureCounts, default (unique assignment), v36
#   6. featureCounts -O --fraction, v36              (shares ambiguous reads)
#
# ATRIP (the gene that contains TREX1) is the internal control; the number of
# read pairs HTSeq puts in __ambiguous is recorded for each run.
set -euo pipefail

W=/media/nighthawk/_mnt_16TB/webtool/validation
IDX=/media/nighthawk/_mnt_3TB/ref_genomes/GRCh38_STAR-2.7.10b_GTF101_oh100
V36=$W/ref/gencode.v36.annotation.gtf
V22=$W/ref/gencode.v22.annotation.gtf
THREADS=48
HT_PROC=16

cd "$W"
mkdir -p bam counts logs

for S in HepG2 sigmoid_colon; do
    R1=fastq/${S}_r1.fastq.gz
    R2=fastq/${S}_r2.fastq.gz
    BAM=bam/${S}_Aligned.out.bam

    if [ ! -s "$BAM" ]; then
        echo "[$(date +%H:%M)] STAR $S"
        STAR --genomeDir "$IDX" \
             --sjdbGTFfile "$V36" --sjdbOverhang 100 \
             --readFilesIn "$R1" "$R2" --readFilesCommand zcat \
             --runThreadN $THREADS \
             --outSAMtype BAM Unsorted \
             --outSAMattributes NH HI AS NM MD \
             --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
             --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
             --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
             --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
             --quantMode GeneCounts \
             --limitSjdbInsertNsj 2000000 \
             --outFileNamePrefix bam/${S}_ > logs/star_${S}.log 2>&1
    fi

    echo "[$(date +%H:%M)] htseq intersection-nonempty v36 $S"
    [ -s counts/${S}_htseq_ine_v36.tsv ] || htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id \
        -m intersection-nonempty -n $HT_PROC "$BAM" "$V36" > counts/${S}_htseq_ine_v36.tsv 2> logs/htseq_ine_v36_${S}.log

    echo "[$(date +%H:%M)] htseq intersection-nonempty v22 $S"
    [ -s counts/${S}_htseq_ine_v22.tsv ] || htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id \
        -m intersection-nonempty -n $HT_PROC "$BAM" "$V22" > counts/${S}_htseq_ine_v22.tsv 2> logs/htseq_ine_v22_${S}.log

    echo "[$(date +%H:%M)] htseq union v36 $S"
    [ -s counts/${S}_htseq_union_v36.tsv ] || htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id \
        -m union -n $HT_PROC "$BAM" "$V36" > counts/${S}_htseq_union_v36.tsv 2> logs/htseq_union_v36_${S}.log

    echo "[$(date +%H:%M)] featureCounts default v36 $S"
    [ -s counts/${S}_fc_default_v36.tsv ] || featureCounts -p --countReadPairs -T $THREADS -t exon -g gene_id \
        -a "$V36" -o counts/${S}_fc_default_v36.tsv "$BAM" > logs/fc_default_${S}.log 2>&1

    echo "[$(date +%H:%M)] featureCounts -O --fraction v36 $S"
    [ -s counts/${S}_fc_fraction_v36.tsv ] || featureCounts -p --countReadPairs -O --fraction -T $THREADS -t exon -g gene_id \
        -a "$V36" -o counts/${S}_fc_fraction_v36.tsv "$BAM" > logs/fc_fraction_${S}.log 2>&1
done

echo "[$(date +%H:%M)] all done"
