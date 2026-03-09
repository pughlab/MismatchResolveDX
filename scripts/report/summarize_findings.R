### summarize_findings.R ###########################################################################
# Summarize findings for key genes.

### FUNCTIONS ######################################################################################
# function to generate a standardized filename
generate.filename <- function(project.stem, file.core, extension, include.date = TRUE) {

	# build up the filename
	file.name <- paste(project.stem, file.core, sep = '_');
	file.name <- paste(file.name, extension, sep = '.');

	if (include.date) {
		file.name <- paste(Sys.Date(), file.name, sep = '_');
		}

	return(file.name);
	}

# function to write session profile to file
save.session.profile <- function(file.name) {

	# open the file
	sink(file = file.name, split = FALSE);

	# write memory usage to file
	cat('### MEMORY USAGE ###############################################################');
	print(proc.time());

	# write sessionInfo to file
	cat("\n### SESSION INFO ###############################################################");
	print(sessionInfo());

	# close the file
	sink();

	}

### PREPARE SESSION ################################################################################
# import command line arguments
library(argparse);

parser <- ArgumentParser();

parser$add_argument('-p', '--project', type = 'character', help = 'PROJECT name');
parser$add_argument('-o', '--output_dir', type = 'character', help = 'path to output directory');
parser$add_argument('-y', '--sample_yaml', type = 'character', help = 'path to sample (BAM) yaml');
parser$add_argument('-c', '--cnas', type = 'character', help = 'path to CNA calls');
parser$add_argument('-g', '--germline', type = 'character', help = 'path to germline MAF file');
parser$add_argument('-s', '--somatic', type = 'character', help = 'path to ensemble (somatic) MAF file');
parser$add_argument('-m', '--methylation', type = 'character', help = 'path to methylation data');
parser$add_argument('-i', '--ichor', type = 'character', help = 'path to ichorCNA TF estimates');

arguments <- parser$parse_args();

###
#arguments$project <- 'UHNBB';
#arguments$sample_yaml <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/UHNBB/DynaCare/DNASeq/GATK/gatk_bam_config_2026-02-09_14-13-48.yaml'
#arguments$output_dir <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/Patient_Reports/SUMMARY';
#arguments$ichor <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/Patient_Reports/data/sWGS/CNAs/2026-02-09_UHNBB_ichorCNA_estimates.tsv'
#arguments$germline <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/Patient_Reports/data/DNASeq/Germline/2026-02-10_UHNBB_mutations_for_cbioportal.tsv';
#arguments$somatic <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/Patient_Reports/data/DNASeq/ENSEMBLE/UHNBB_ensemble_mutation_data.tsv';
#arguments$cnas <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/Patient_Reports/data/DNASeq/CNAs/2026-02-10_UHNBB_panelCN.mops_output.tsv';
#arguments$methylation <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/Patient_Reports/data/EMSeq/Methylation/2025-10-10_MultiMMR_methylation_per_target_region.tsv'
###


# load libraries
library(yaml);

# what's the date?
date <- Sys.Date();

# create (if necessary) and move to output directory
output.dir <- arguments$output_dir;
if (!dir.exists(output.dir)) {
	dir.create(output.dir);
	}

setwd(output.dir);

# indicate key genes
classic.genes <- unlist(strsplit('MLH1|MSH2|MSH6|PMS2','\\|'));
ancillary.genes <- unlist(strsplit('EPCAM|MLH3|POLD1|POLE|MUTYH|TP53|BRAF','\\|'));

genes.of.interest <- c(classic.genes, ancillary.genes);

### READ DATA ######################################################################################
# get data
if (!is.null(arguments$ichor)) {
	tf.data <- read.delim(arguments$ichor);
	} else {
	tf.data <- NULL;
	}

if (!is.null(arguments$germline)) {
	germline.data <- read.delim(arguments$germline, comment.char = '#');
	} else {
	germline.data <- NULL;
	}

if (!is.null(arguments$somatic)) {
	somatic.data <- read.delim(arguments$somatic, comment.char = '#');
	} else {
	somatic.data <- NULL;
	}

if (!is.null(arguments$cnas)) {
	cna.data <- read.delim(arguments$cnas, stringsAsFactors = FALSE);
	} else {
	cna.data <- NULL;
	}

if (!is.null(arguments$methylation)) {
	methylation.data <- read.delim(arguments$methylation);
	} else { 
	methylation.data <- NULL;
	}

# parse sample information from yaml file
project.yaml <- read_yaml(arguments$sample_yaml);

sample.info <- as.data.frame(matrix(nrow = 0, ncol = 3));
colnames(sample.info) <- c('Patient','Sample','Type');

for (patient in names(project.yaml)) {
	normals <- names(project.yaml[[patient]]$normal);
	tumours <- names(project.yaml[[patient]]$tumour);
	sample.info <- rbind(sample.info, 
		data.frame(Patient = rep(patient, length(normals)), Sample = normals, Type = rep('normal', length(normals))),
		data.frame(Patient = rep(patient, length(tumours)), Sample = tumours, Type = rep('tumour', length(tumours)))
		);
	}

sample.info$Patient <- as.character(sample.info$Patient);
sample.info$Sample <- as.character(sample.info$Sample);

sample.info <- sample.info[order(sample.info$Patient, sample.info$Sample),];


### FORMAT DATA ####################################################################################
# format somatic SNV/INDEL data
if (!is.null(somatic.data)) {
	somatic.data <- somatic.data[which(somatic.data$Hugo_Symbol %in% genes.of.interest),];
	somatic.data$VAF <- somatic.data$t_alt_count / somatic.data$t_depth;
	somatic.data$ClinVar <- sapply(somatic.data$CLIN_SIG, function(i) {
		if (is.na(i)) { return(NA) 
			} else {
			tmp <- unlist(strsplit(i,','));
			x <- if (any(tmp == 'pathogenic')) { 'pathogenic'
				} else if (any(tmp == 'likely_pathogenic')) { 'likely_pathogenic'
				} else if (any(tmp == 'likely_benign')) { 'likely_benign'
				} else if (any(tmp == 'benign')) { 'benign'
				} else if (any(tmp %in% c('uncertain_significance','conflicting_interpretations_of_pathogenicity'))) { 'VUS'
				} else { NA }
			return(x);
			}
		});

	keep.fields <- c('Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification','HGVSc','HGVSp_Short','ClinVar','VAF');

	somatic.data <- somatic.data[!grepl('benign', somatic.data$ClinVar),keep.fields];
	somatic.data <- somatic.data[!grepl('Flank|UTR|Intron|Silent', somatic.data$Variant_Classification),];
	}

if (!is.null(germline.data)) {
	germline.data <- germline.data[which(germline.data$Hugo_Symbol %in% genes.of.interest),];
	germline.data$VAF <- germline.data$t_alt_count / germline.data$t_depth;
	germline.data$ClinVar <- sapply(germline.data$CLIN_SIG, function(i) {
		if (is.na(i)) { return(NA) 
			} else {
			tmp <- unlist(strsplit(i,','));
			x <- if (any(tmp == 'pathogenic')) { 'pathogenic'
				} else if (any(tmp == 'likely_pathogenic')) { 'likely_pathogenic'
				} else if (any(tmp == 'likely_benign')) { 'likely_benign'
				} else if (any(tmp == 'benign')) { 'benign'
				} else if (any(tmp %in% c('uncertain_significance','conflicting_interpretations_of_pathogenicity'))) { 'VUS'
				} else { NA }
			return(x);
			}
		});

	keep.fields <- c('Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification','HGVSc','HGVSp_Short','ClinVar','VAF');

	germline.data <- germline.data[grepl('pathogenic', germline.data$ClinVar),keep.fields];
	}

# format CNAs
# indicate region of EPCAM-MSH2 deletion
for (smp in unique(cna.data$Sample)) {
	idx1 <- which(cna.data$Sample == smp & cna.data$Gene == 'EPCAM')[1];
	idx2 <- which(cna.data$Sample == smp & cna.data$Gene == 'MSH2')[1];

	cna.data[idx1:idx2,]$Gene <- 'EPCAM'
	}

# exclude 'other' sites (these are generally MSI sites)
cna.data <- cna.data[!grepl('other', cna.data$Exon),];

# extract genes of interest
cna.data <- cna.data[which(cna.data$Gene %in% genes.of.interest),];
cna.data$CN <- as.numeric(gsub('CN','',as.character(cna.data$CN)));
cna.data[which(cna.data$lowQual == 'lowQual'),]$CN <- NA;

# summarize gene-level CN status using weighted mean
gene.cnas <- aggregate(
	CN ~ Sample + Gene,
	cna.data,
	mean,
	na.rm = TRUE
	);

gene.cnas <- gene.cnas[which(gene.cnas$CN > 2.8 | gene.cnas$CN < 1.7),];

# format methylation
methyl.samples <- colnames(methylation.data)[7:ncol(methylation.data)];
methyl.gene.data <- reshape(
	methylation.data,
	direction = 'long',
	varying = list(7:ncol(methylation.data)),
	timevar = 'Sample',
	times = methyl.samples,
	v.names = 'Fraction'
	);

methyl.gene.data$Gene <- sapply(methyl.gene.data$V4, function(i) { unlist(strsplit(i,'_'))[1] } );
methyl.gene.data <- methyl.gene.data[which(methyl.gene.data$Gene %in% genes.of.interest),];
methyl.gene.data$Sample <- gsub('\\.','-',methyl.gene.data$Sample);


# initiate data objects
if (!is.null(arguments$ichor)) {

	master.matrix <- merge(
		sample.info[which(sample.info$Type != 'normal'),1:2],
		tf.data[,c('ID','Tumour.Fraction')],
		by.x = 'Sample',
		by.y = 'ID',
		all.x = TRUE
		);

	} else {

	master.matrix <- sample.info[which(sample.info$Type != 'normal'),1:2];
	master.matrix$Tumour.Fraction <- NA;

	}


# add mutations for each gene
for (gene in genes.of.interest) {

	tmp <- master.matrix;
	tmp[,c('germline','somatic','methylation','copynumber')] <- NA;

	df1 <- unique(germline.data[which(germline.data$Hugo_Symbol == gene),]$Tumor_Sample_Barcode);
	df2 <- unique(somatic.data[which(somatic.data$Hugo_Symbol == gene),]$Tumor_Sample_Barcode);
	df3 <- methyl.gene.data[which(methyl.gene.data$Gene == gene & methyl.gene.data$Fraction > 0.1),]$Sample;
	df4.up <- gene.cnas[which(gene.cnas$Gene == gene & gene.cnas$CN > 2),]$Sample;
	df4.dn <- gene.cnas[which(gene.cnas$Gene == gene & gene.cnas$CN < 2),]$Sample;

	if (length(df1) > 0) {
		tmp[which(tmp$Sample %in% df1),]$germline <- 'germline';
		}
	if (length(df2) > 0) {
		tmp[which(tmp$Sample %in% df2),]$somatic <- 'somatic';
		}
	if (length(df3) > 0) {
		tmp[which(tmp$Sample %in% df3),]$methylation <- 'methylation';
		}
	if (length(df4.up) > 0) {
		tmp[which(tmp$Sample %in% df4.up),]$copynumber <- 'cn.gain';
		}
	if (length(df4.dn) > 0) {
		tmp[which(tmp$Sample %in% df4.dn),]$copynumber <- 'cn.loss';
		}

	tmp[,gene] <- apply(tmp[,c('germline','somatic','methylation','copynumber')],1,function(i) {
		paste(na.omit(i), collapse = '/') });

	master.matrix <- merge(
		master.matrix,
		tmp[,c('Sample',gene)],
		all.x = TRUE
		);
	}

# save results
write.table(
	master.matrix[,c('Patient','Sample','Tumour.Fraction',genes.of.interest)],
	file = generate.filename(arguments$project, '_final_summary','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename(arguments$project, 'SummarizeFindings_SessionProfile','txt'));
