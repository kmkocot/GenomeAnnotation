#https://github.com/Gaius-Augustus/BRAKER
#https://www.biostars.org/p/411101/
#Make sure to remove line breaks in sequences in protein_evidence.fas
#Make sure to remove whitespace in genome fasta header


#Set correct paths for GeneMark and Augustus 3.4.0
export AUGUSTUS_BIN_PATH=/home/wirenia/bin/Augustus-3.5.0/bin
export AUGUSTUS_SCRIPTS_PATH=/home/wirenia/bin/Augustus-3.5.0/scripts
export AUGUSTUS_CONFIG_PATH=/home/wirenia/bin/Augustus-3.5.0/config

#Define variables
GENOME=out_JBAT.FINAL.chr.fa
FORWARD_READS=fwd.fq
REVERSE_READS=rev.fq
PROTEIN_EVIDENCE=protein_evidence.fas
AUGUSTUS_SPECIES=`date +%Y-%m-%d_%H.%M_$GENOME`
CORES=30


#Quality filter and trim transcriptome reads
trim_galore --cores $CORES --fastqc --quality 30 --length 50 --adapter file:/home/wirenia/databases/adapters.fasta --paired $FORWARD_READS $REVERSE_READS
gunzip *_val_1.fq.gz
gunzip *_val_2.fq.gz


#Prepare protein_evidence.fas file - expects protein fasta files with a .pep extension to be in a directory called "protein_evidence" 
cat ./protein_evidence/*.pep > ./protein_evidence.fas
awk '!/^>/ { printf "%s", $0; n = "\n" } 
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' protein_evidence.fas > protein_evidence.fas.nent
rm protein_evidence.fas
rename 's/.nent//g' *.nent
sed -i 's/ .\+//g' protein_evidence.fas #Cleans up fasta headers

#Make RepeatModeler database
BuildDatabase -engine rmblast -name repeats $GENOME


#Run RepeatModeler - use only ~16 cores max for large (e.g., human-sized or bigger) genomes
RepeatModeler -pa $CORES -engine rmblast -LTRStruct -database repeats 2>&1 | tee repeatmodeler.log


#Run RepeatMasker - use only ~16 cores max for large (e.g., human-sized or bigger) genomes
cp ./*/consensi.fa.classified .
RepeatMasker -parallel $CORES -engine rmblast -gccalc -lib ./consensi.fa.classified $GENOME -xsmall


#Make hisat2 index for masked genome
hisat2-build $GENOME".masked" index


#Map transcriptome reads to masked genome with hisat2
hisat2 -x index -1 *_val_1.fq -2 *_val_2.fq -S Aligned.out.sam


#Make BAM from SAM and get rid of SAM and hisat2 index files:
samtools view -bS Aligned.out.sam > RNAseq.bam
rm -rf *.sam
rm -rf *.ht2


#Run BRAKER
#Don't use --crf
#perl /home/wirenia/bin/BRAKER-3.0.8/scripts/braker.pl --threads $CORES --softmasking --makehub --email kmkocot@ua.edu --gff3 --species $AUGUSTUS_SPECIES --prot_seq protein_evidence.fas --bam RNAseq.bam --genome $GENOME".masked"
#perl /home/wirenia/bin/BRAKER-3.0.8/scripts/braker.pl --threads $CORES --softmasking --makehub --email kmkocot@ua.edu --gff3 --species $AUGUSTUS_SPECIES --prot_seq protein_evidence.fas --genome $GENOME".masked"
perl /home/wirenia/.local/bin/BRAKER-3.0.8/scripts/braker.pl --threads $CORES --softmasking --makehub --email kmkocot@ua.edu --gff3 --prot_seq=protein_evidence.fas --bam=RNAseq.bam --species $AUGUSTUS_SPECIES --genome $GENOME".masked"


#Run BUSCO
cd braker
busco --cpu $CORES -m proteins -l metazoa_odb10 --download_path /home/wirenia/databases -i braker.aa -o BUSCO_braker.aa
busco --cpu $CORES -m proteins -l metazoa_odb10 --download_path /home/wirenia/databases -i augustus.hints.aa -o BUSCO_augustus.hints.aa

