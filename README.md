# PEPG pipeline's description was listed in the PEPG.pdf file in destails. 
The first step is the MAKER2-Legacy-Annotation to add AED Quality metrics to the prior Annotation.

# 	Description

With all new RNA and Protein homology evidence, MAKER2 offers a way to accomplish measuring the quality of gene prediction from any preexisting annotation by adding the same quality metrics, like Annotation Edit Distance (AED) and mRNA Quality index (QI). AED is a number between 0 and 1, with an AED score of 0 indicating a perfect prediction with available evidence and a value of 1 without any evidence supporting the annotated gene model. A prior annotation file in GFF3 or GTF format is the start point here at first.

# (Optional, skip this step if you have it ready) And then, in order to feed the PEPG pipeline, the RNA transcriptome evidence are prepared here: 
The reference genome based transcriptome assemble pipeline of RNA-seq reads can be found in https://github.com/didiren/transcriptome_assembly.git

# (Optional, skip this step if you have it ready) This is a evaluation pipeline by comparison, so if you don't have a new.annotation.gff at first, then start from here:
MAKER2-two-pass genome annotation pipeline's description was listed in the MAKER2-two-pass pipeline.pdf file in destails if you want to re-annotate your genome with the RNA-seq evidence and updated protein evidence. 

# 1.1	Running MAKER2-Legacy 

With the same strategy like normal MAKER2 one-pass, we need to edit the three control files for MAKER2 to get access to the genome, RNA, Proteins evidence. However, there is no need to utilize any gene predictor because preexisting annotation provided gene models for MAKER2 to exam with the evidence it has.
#Details: PEPG.pdf


# 1.2 Running InterproScan

The scond step is to run the InterproScan append the Interpro conserved domains ID to each predicted gene from the prior annotation, so the comparison can be in both functional and strutural of each gene between the prior annotation and the newer version.
#Details: PEPG.pdf

# 1.3 Running PEPG.py
The last step is to run the PEPG.py, which is designed  to compare the prior annotation and the newer annotation in the angles of the accuracy, gene's structure and function. It requires you to provides the two GFF files from the two annotations.

#Commandline: python PEPG.py -g1 prior.annotation.gff -g2 new.annotation.gff




