
#!/bin/bash

#
# eppi - analysis pipeline for viral metagenomics
#

PATH=.:$PATH
source ~/ma/lab/15-eppi/eppi.config
reads_file=$1

mkdir out

#fastQC
mkdir out/fastqc
mkdir out/fastqc/before-filtering
fastqc --threads $threads --outdir out/fastqc/before-filtering $reads_file

# filter reads
fastq_quality_filter -Q33 -q 20 -p 90 -i $reads_file -o out/qc_filtered.fastq
reads_file="out/qc_filtered.fastq"

#fastQC
mkdir out/fastqc/after-filtering
fastqc --threads $threads --outdir out/fastqc/after-filtering $reads_file

# Vrap assembly
mkdir out/vrap
vrap.py --noblast -g -1 $reads_file -t $threads -o out/vrap
cp out/vrap/vrap_contig.fasta out/contigs.fasta

# assembly check
mkdir out/metaquast
~/ma/scr/quast-4.5/metaquast.py --threads $threads --output-dir out/metaquast out/contigs.fasta

# Mapping reads against contigs
mkdir out/read-cov
hisat2-build -p $threads out/contigs.fasta out/read-cov/contigs # build index with contigs
hisat2 -p $threads -x out/read-cov/contigs -U $reads_file -S out/read-cov/reads-against-contigs.sam
sam2bam out/read-cov/reads-against-contigs.sam
rm out/read-cov/reads-against-contigs.sam

# read coverage for contigs
coverage-density.py out/read-cov/reads-against-contigs.sam.bam out/read-cov/

# Nucleotide read alignment
mkdir out/centrifuge
centrifuge -p $threads --time --met-file out/centrifuge/metrics.txt -S out/centrifuge/hits.tsv \
	--report-file out/centrifuge/centrifuge_report.tsv --min-hitlen $centrifuge_min_hit_length \
	-x $centrifuge_index -U $reads_file
cat out/centrifuge/hits.tsv | tail -n +2 | cut -f 3 | sort | uniq -c | sed -E s/^" "+//g | sed s/" "/","/g > out/centrifuge/taxon-counts.csv
centrifuge_tax_count_to_lineage.py out/centrifuge/taxon-counts.csv $tax_db out/centrifuge/lineage-counts.csv


# Protein read alignment
mkdir out/diamond
diamond blastx -p $threads --db $diamond_db --query $reads_file --top 5 --daa out/diamond/blastx-reads-hits
diamond view --daa out/diamond/blastx-reads-hits.daa --outfmt 6 qseqid sseqid sstart send evalue btop --out out/diamond/blastx-reads-hits.tsv


# Nucleotide contig alignment
mkdir out/blast
blastn -db $blast_contig_nt_db -query out/contigs.fasta -outfmt "6 qseqid sseqid sstart send evalue btop" -out out/blast/blastn-contigs-nt.tsv  -num_threads $threads
lca-blast.py out/blast/blastn-contigs-nt.tsv $tax_db

# Protein contig alignment
diamond blastx -p $threads --db $diamond_db --query out/contigs.fasta --daa out/diamond/blastx-contigs-hits
diamond view --daa out/diamond/blastx-contigs-hits.daa --outfmt 6 qseqid sseqid sstart send evalue btop --out out/diamond/blastx-contigs-hits.tsv
lca-blast.py out/diamond/blastx-contigs-hits.tsv $tax_db


# HMM contig alignment
mkdir out/hmm
contigs_file="out/contigs.fasta"
for frame in {-1,-2,-3,1,2,3}; do
	orf_seq=out/hmm/contigs.$frame.fna
	transeq -sequence $contigs_file -frame $frame -outseq $orf_seq
	hmmsearch --domtblout $orf_seq-dom-hits.tsv --cpu $threads $hmm_db $orf_seq
done

# HMM dom out to species taxon hit table
for hmm_dom_out in $(ls out/hmm/*dom-hits*); do
	hmm-dom-out-to-tax-hit-table.py $hmm_dom_out $hmm_dom_out.taxhits.csv $tax_db
done

# HMM lca
lca-hmm.py $tax_db "out/hmm/*.taxhits.csv" out/hmm/lca_df.csv

# Mash
mkdir out/mash
mash screen -p $threads -i 0.8 -v 0.05 $mash_db $reads_file > out/mash/reads-mash-screen.tsv

# deep learning
mkdir out/deep-learning
classify-dnn.py out/contigs.fasta $deep_learning_theta > out/deep-learning/dnn-hits.tsv

# combine contig annotations
annotate-contigs.py

# general stats
#stats.sh
