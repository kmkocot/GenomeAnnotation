These scripts run all the steps needed to:
Quality trim and filter transcriptome reads with trim_galore
Identify repeats in your genome de novo with RepeatModeler
Mask repeats with RepeatMasker
Run ProtHint to generate hints from translated transcriptomes (or do this as part of the BRAKER3 pipeline)
Map trimmed/filtered reads to the genome with STAR/HISAT2
Run BRAKER to annotate the genome

Some relevant reading:
https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
http://www.repeatmasker.org/RepeatModeler/
http://www.repeatmasker.org/RepeatMasker/
https://www.biostars.org/p/411101/
https://github.com/gatech-genemark/ProtHint
https://github.com/Gaius-Augustus/BRAKER/

Note that these scripts are written for Kevin Kocot's workstation and some paths are hard-coded (search for "wirenia" to find them)
