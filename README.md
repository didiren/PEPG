# PEPG pipeline's description was listed in the PEPG.pdf file in destails. 
The first step is the MAKER2-Legacy-Annotation to add AED Quality metrics to the prior Annotation.

1.1	Description

With all new RNA and Protein homology evidence, MAKER2 offers a way to accomplish measuring the quality of gene prediction from any preexisting annotation by adding the same quality metrics, like Annotation Edit Distance (AED) and mRNA Quality index (QI). AED is a number between 0 and 1, with an AED score of 0 indicating a perfect prediction with available evidence and a value of 1 without any evidence supporting the annotated gene model. A preexisting annotation file in GFF3 format and the same kinds of evidence we used above for consistency were used here.

1.2	Running MAKER2-Legacy 

With the same strategy like normal MAKER2 one-pass, we need to edit the three control files for MAKER2 to get access to the genome, RNA, Proteins evidence. However, there is no need to utilize any gene predictor because preexisting annotation provided gene models for MAKER2 to exam with the evidence it has.

The scond step is to run the InterproScan append the Interpro conserved domains ID to each predicted gene from the prior annotation, so the comparison can be in both functional and strutural of each gene between the prior annotation and the newer version.


The last step is to run the PEPG.py, which is designed  to compare the prior annotation and the newer annotation in the angles of the accuracy, gene's structure and function. It requires you to provides the two GFF files from the two annotations.
