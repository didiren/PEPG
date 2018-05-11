# PEPG
The first step is the MAKER2-Legacy-Annotation to add AED Quality metrics to the prior Annotation


1.1	Description

With all new RNA and Protein homology evidence, MAKER2 offers a way to accomplish measuring the quality of gene prediction from any preexisting annotation by adding the same quality metrics, like Annotation Edit Distance (AED) and mRNA Quality index (QI). AED is a number between 0 and 1, with an AED score of 0 indicating a perfect prediction with available evidence and a value of 1 without any evidence supporting the annotated gene model. A preexisting annotation file in GFF3 format and the same kinds of evidence we used above for consistency were used here.

1.2	Running MAKER2-Legacy 

With the same strategy like normal MAKER2 one-pass, we need to edit the three control files for MAKER2 to get access to the genome, RNA, Proteins evidence. However, there is no need to utilize any gene predictor because preexisting annotation provided gene models for MAKER2 to exam with the evidence it has.

1.2.1	Editing three control files
With three control files existing in didiren/ folder, the maker_opts.ctl file was edited to specify all inputs, like Genome sequence, RNA evidence, Protein evidence and the preexisting GFF3 file.

genome=/work/didiren/Cparasitica.genome.fasta
organism_type=eukaryotic
est=/work/didiren/assembled.transcriptome.fasta
protein=/work/didiren/uniprot_sprot.fasta.gz
model_gff=/work/didi/Cparasiticav2.GeneCatalog20091217.gff3
est2genome=0
protein2genome=0


The option description here:
#genome---to give genome sequences
#organism_type           
to indicate the type of organism
#est                
to give RNA evidence sequences 
#protein                    
to give protein evidence sequences 
#model_gff   
to give external preexisting annotated gene models 
#est2genome          
to infer gene predictions directly from RNA sequence, 1=yes,     0=no
#protein2genome     
    to infer gene predictions directly from protein homology, 1=yes, 0=no


1.2.2	To run MAKER2-Legacy 

Once maker_opts.ctl have been edited, run Maker2 by creating this bash script file maker_legacy.sh in the /work/didiren folder and submit it to the sever.

#!/bin/bash
#PBS -N Run_Maker_FirstPass
#PBS -l nodes=1:ppn=20
#PBS -l walltime=48:00:00

#========================= SETUP MAKER =========================
swsetup () { eval `/usr/erc-share/etc/swsetup/swsetup.pl $*`; }
swsetup maker
swsetup augustus
PERL5LIB=/usr/local/igbb/maker/lib/perl5:$PERL5LIB
PERL5LIB=/usr/local/igbb/genemark-es-et_4.30/lib/perl5:$PERL5LIB
#---------------------------------------------------------------

# CHANGE DIRECTORY TO WHERE JOB WAS SUBMITTED
cd /work/didiren

# RUN MAKER2 WITH MAX NUMBER OF THREADS
maker -c $PBS_NUM_PPN -base maker2_legacy
# -c        to use multiple processors
# -base     to use Set the base name MAKER uses to save output files



Now submit the job to the sever by typing this command at the directory of where your maker_legacy.sh file is.

qsub maker_legacy.sh

1.2.3	MAKER2-Legacy outputs

Now in the current working directory, MAKER2 have created a folder named maker2_legacy.maker.output. Inside this directory, maker2_legacy_master_datastore_index.log need to be looked into first to make sure that all scaffolds were taken successfully. As we have done before, they were combined to generate an integrated gff, proteins.fasta, transcripts.fasta first by creating this maker_combine1.sh bash scripts file.

#!/bin/bash
#PBS -N Run_MAKER_COMBINE1
#PBS -l nodes=1:ppn=20 
#PBS -l walltime=48:00:00 

#========================= SETUP MAKER =========================
swsetup () { eval `/usr/erc-share/etc/swsetup/swsetup.pl $*`; }
swsetup maker
swsetup augustus
PERL5LIB=/usr/local/igbb/maker/lib/perl5:$PERL5LIB
PERL5LIB=/usr/local/igbb/genemark-es-et_4.30/lib/perl5:$PERL5LIB
#---------------------------------------------------------------

# CHANGE DIRECTORY TO WHERE JOB WAS SUBMITTED
cd /work/didiren/

base=maker2_legacy
# to merge all gff3 frile from each scaffold to one gff3 file
gff3_merge -d $base.maker.output/$base\_master_datastore_index.log
#to merge all fasta frile from each scaffold to one fasta file
fasta_merge -d $base.maker.output/$base\_master_datastore_index.log

Now submit the job to the sever by typing this command at the directory of where your maker_combine1.sh file is.

qsub maker_combine1.sh

2	AED Comparison between preexisting and new annotation
2.1	Extraction of AED score from the GFF3 file. 
In a GFF3 format file, there are 9 columns, containing seqid, source, type, start, end, score, strand, phase and attributes, separately. In the ninth column, gene ID, gene name, AED score, etc were presented by a semicolon-separated list. A R script were written to extract AED score using regulator expression function as followed.


 
