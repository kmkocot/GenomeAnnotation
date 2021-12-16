#!/bin/bash

#This script runs all the steps needed to:
#1) Identify repeats in your genome de novo with RepeatModeler
#2) Mask repeats with RepeatMasker
#3) Quality trim and filter transcriptome reads with trim_galore
#4) Map trimmed/filtered reads to the genome with STAR
#5) Run ProtHint to generate hints from translated transcriptomes
#5) Run braker.pl to annotate the genome

#Relevant manuals and reading:
#https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
#http://www.repeatmasker.org/RepeatModeler/
#http://www.repeatmasker.org/RepeatMasker/
#https://www.biostars.org/p/411101/ - create a de novo repeat library
#https://github.com/gatech-genemark/ProtHint
#https://github.com/Gaius-Augustus/BRAKER/

#Specify input file names, number of cores to use, etc:
EMAIL=kmkocot@ua.edu #Needed if you want to make a UCSC genome browser with braker.pl
CORES=38 #Number of CPU cores to use; don't use too many for trim_galore
SPECIES=Hanleya_hanleyi #Species name (don't use spaces)
GENOME=final.purged.fa.PolcaCorrected.fa #Genome assembly filename
RAW_TR_1=Hanleya_nagelfar_mantle_1.fastq #Left transcriptome reads
RAW_TR_2=Hanleya_nagelfar_mantle_2.fastq #Right transcriptome reads

#Leave these alone
TRIM_TR_1=`echo $RAW_TR_1 | awk -F "." '{print $1}'`
TRIM_TR_2=`echo $RAW_TR_2 | awk -F "." '{print $1}'`
export HERE=`pwd`

#Trim and filter transcriptome data:
#Don't use more than 8 cores.
gunzip Hanleya_nagelfar_mantle_1.fastq.gz
gunzip Hanleya_nagelfar_mantle_2.fastq.gz
/home/scutopus/bin/TrimGalore-0.6.7/trim_galore --cores 2 -q 30 --illumina --length 50 --trim-n -o $HERE --paired $RAW_TR_1 $RAW_TR_2

#Make RepeatModeler database:
BuildDatabase -engine rmblast -name $SPECIES"_repeats" $GENOME

#Run RepeatModeler:
RepeatModeler -pa $CORES -engine rmblast -LTRStruct -database $SPECIES"_repeats" 2>&1 | tee repeatmodeler.log

#Run RepeatMasker:
#GC content can be optimized for your organism if you specify the -gc flag
RepeatMasker -parallel $CORES -engine rmblast -lib ./*/consensi.fa.classified Hanleya_hanleyi_repeats $GENOME -xsmall

#Run ProtHint:
#Make sure there are no line breaks in the evidence file!
prothint.py --threads $CORES --cleanup --evalue 1e-25 $GENOME".masked" protein_evidence.fas

#Generate STAR genome:
#Make sure to use the repeat-masked genome.
mkdir star_genome
STAR --runThreadN $CORES --runMode genomeGenerate --genomeDir star_genome --genomeFastaFiles $GENOME".masked"

#Run STAR mapping:
#Specify the ram limit for your machine.
STAR --runThreadN $CORES --limitGenomeGenerateRAM 63000000000 --chimSegmentMin 50 --outFilterType BySJout --genomeDir star_genome --readFilesIn $TRIM_TR_1"_val_1.fq" $TRIM_TR_2"_val_2.fq"

#Make BAM from SAM and get rid of SAM:
samtools view -bS Aligned.out.sam > RNAseq.bam
rm -rf *.sam

#Run BRAKER2:
#Note that --UTR=on may not be compatible with adding hints from ProtHint(?). You can always add UTRs later with a re-run of braker.pl.
#You must use Augustus 3.3.3 (not 3.4.0) BUT replace the Augustus 3.3.3 joingenes with the one from 3.4.0.
#Don't forget to update the --species flag if you are training Augustus and did a previous run with the same species name.
braker.pl --cores $CORES --etpmode --verbosity=4 --softmasking --filterOutShort --crf --makehub --email $EMAIL --gff3 --species $SPECIES --genome $GENOME".masked" --bam RNAseq.bam --hints=prothint_augustus.gff

#Assess the output:
cd braker
echo "Number of gene models in final prediction:"
grep -c ">" augustus.hints.aa
echo
busco -c 20 -m protein --long -l metazoa_odb10 -i augustus.hints.aa -o BUSCO_augustus.hints.aa
busco -c 20 -m transcriptome --long -l metazoa_odb10 -i augustus.hints.codingseq -o BUSCO_augustus.hints.codingseq



