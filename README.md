# MismatchResolveDx (version 1.0.0)

## Introduction
This is a collection of pipelines to be used for processing the MismatchResolveDx (formarly MultiMMR) target panel (including sWGS, targeted DNA-seq and targeted EM-Seq). These pipelines are based on the PughLab Pipeline-Suite (https://github.com/pughlab/pipeline-suite.git).
- pughlab_multimmr_pipeline.pl: master pipeline; submits calls to pughlab_dnaseq_pipeline.pl and pughlab_emseq_pipeline.pl
- pughlab_dnaseq_pipeline.pl: performs adapter trimming, alignment, qc and variant-calling on sWGS and targeted DNA-seq
- pughlab_emseq_pipeline.pl: performs adapter trimming, alignment, qc and methylation calling on targeted EM-Seq data


Start by creating a clone of the repository:

<pre><code>cd /path/to/some/directory
git clone git@github.com:pughlab/MismatchResolveDX.git
</code></pre>

## Dependencies
Pipeline-Suite was developed using perl and is currently compatible using the Slurm workload manager.

In addition, all of the desired tools (ie, BWA, GATK, etc) must be installed and be accessible using the 'module load \<tool\/version\>' syntax.

Tools required for processing/analysis of sWGS (all versions shown are those used during pipeline development): 
- perl (v5.40.0) and R (v4.1.0)
- samtools (v1.20) and picardtools (v2.10.9)
- fastp (v0.23.1) for adapter trimming
- fastqc (v0.11.5) for QC of trimmed fastq files
- BWA (v0.7.15) for alignment, with picardtools for mark duplicates
- picardtools and GATK (v4.6.0.0) to collect BAM QC metrics
- hmmcopy_utils readCounter (v170718) to generate WIG files for ichorCNA
- R packages required for ichorCNA including ichorCNA (v0.3.2), HMMcopy (v1.34.0)

Tools required for processing/analysis of targeted DNA-Seq (all versions shown are those used during pipeline development): 
- perl (v5.40.0) and R (v4.1.0)
- samtools (v1.20), picardtools (v2.10.9), vcftools (v0.1.15)
- trim_galore (v0.6.6) for adapter trimming
- fastqc (v0.11.5) for QC of trimmed fastq files
- BWA (v0.7.15)for alignment, with picardtools for mark duplicates
- GATK (v3.8) for indel realignment and base quality score recalibration
- picardtools (v2.10.9) and GATK (v4.6.0.0) to collect BAM QC metrics
- HaplotypeCaller/GenotypeGVCFs (GATK v3.8) for germline variant calling
- CPSR (v0.6.1) and PCGR (v0.9.1), VEP (v98) and vcf2maf (v1.6.17) for variant annotations
- MuTect2 (GATK v3.8), Pindel (v0.2.5b8), SomaticSniper (v1.0.5.0)/bam-readcount (v0.7.4) and Vardict (java;v1.7.0) for somatic variant calling
- Delly (v0.8.1), novoBreak (v1.1) and Mavis (v2.2.5) for SV detection
- R packages required for panelCN.mops including panelcn.mops (v1.14.0), CopyNumberPlots (v) and CopyNumberPlots (v1.8.0)

Tools required for processing/analysis of targeted EM-Seq (all versions shown are those used during pipeline development): 
- perl (v5.40.0), R (v4.1.0), python3 (v3.10.9)
- samtools (v1.20) and picardtools (v2.10.9)
- trim_galore (v0.6.6) for adapter trimming
- fastqc (v0.11.5) for QC of trimmed fastq files
- BWA (v0.7.15), BWA-METH (v0.2.7) for alignment, with sambamba (v0.7.0) for mark duplicates
- picardtools (v2.10.9) to collect BAM QC metrics
- MethylDackel (v0.6.1) to collect per-site methylation estimates
- R packages required to format methylation data including BiocParallel (v1.28.3) and bsseq (v1.28.0)

Other R packages including:
- optparse, argparse, plyr, yaml
- GenomicRanges (v1.46.1) and GenomeInfoDb (v1.30.1)
- org.Hs.eg.db (v3.13.0) and AnnotationDbi (v1.54.0)
- BSgenome.Hsapiens.UCSC.hg38 (v1.4.3) and/or BSgenome.Hsapiens.UCSC.hg19
- TxDb.Hsapiens.UCSC.hg38.knownGene (v3.13.0) and/or TxDb.Hsapiens.UCSC.hg19.knownGene
- for visualizations: BPG (https://CRAN.R-project.org/package=BoutrosLab.plotting.general), UpSetR (https://cran.r-project.org/web/packages/UpSetR/index.html) and RCircos (https://github.com/hzhanghenry/RCircos)

## Workflow
See the [wiki](https://github.com/pughlab/pipeline-suite/wiki) for detailed instructions on how to set up a new project and run the pipeline.
