# TetV RNA

This serves as a repository of codes to analyze Tetraselmis sp - TetV1 metatranscriptome. 

**VIRUS ANALYSIS**

1. **Trimmomatic and FastQC**  
2. **Tuxedo protocol** (Genome-guided transcriptome assembly)
      a. **HISAT2**. RNA-seq PE read alignment.
      b. **StringTie**. Assembly, quantification, merge all transcripts.
      c. **Ballgown.** Calculate FPKM.

**HOST ANALYSIS**

1. **Bowtie2** (remove virus, TetV-1, reads)
2. **Trinity** (de novo assembly)
3. **BUSCO** (Chlorophyta single copy genes)
4. **TransRate** (filter the bad assemblies) 
5. **RSEM** (Gene/Transcript quantitation)
6. **DESeq2** in Trinity (Differential Expression)
7. **TransDecoder** (coding sequence identification)
8. **Trinotate** (functional annotation)
9. **topGO** (gene ontology enrichment)
