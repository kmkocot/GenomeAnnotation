These scripts run all the steps needed to:
1) Quality trim and filter transcriptome reads with trim_galore,
2) Identify repeats in your genome de novo with RepeatModeler,
3) Mask repeats with RepeatMasker,
4) Run ProtHint to generate hints from translated transcriptomes (or do this as part of the BRAKER3 pipeline),
5) Map trimmed/filtered reads to the genome with STAR/HISAT2,
6) Run BRAKER to annotate the genome


Some relevant reading:
https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
http://www.repeatmasker.org/RepeatModeler/
http://www.repeatmasker.org/RepeatMasker/
https://www.biostars.org/p/411101/
https://github.com/gatech-genemark/ProtHint
https://github.com/Gaius-Augustus/BRAKER/
https://github.com/Gaius-Augustus/TSEBRA

Note that these scripts are written for Kevin Kocot's workstation and some paths are hard-coded (search for "wirenia" to find them)
