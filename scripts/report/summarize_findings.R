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
parser$add_argument('-i', '--ichor', type = 'character', help = 'path to ichorCNA TF estimates');

arguments <- parser$parse_args();

###
#arguments$project <- 'UHNBB';
#arguments$sample_yaml <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/UHNBB/DynaCare/DNASeq/GATK/gatk_bam_config_2026-02-09_14-13-48.yaml'
#arguments$output_dir <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/Patient_Reports/SUMMARY';
#arguments$ichor <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/Patient_Reports/data/sWGS/CNAs/2026-02-09_UHNBB_ichorCNA_estimates.tsv'
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

# find and read in mutation data
mutation.files <- list.files(path = '../Mutations', pattern = arguments$project, full.names = TRUE);
germline.files <- rev(sort(mutation.files[grepl('germline_mutations', mutation.files)]));
somatic.files <- rev(sort(mutation.files[grepl('somatic_mutations', mutation.files)]));

if (length(germline.files) > 0) {
	germline.data <- read.delim(germline.files[1]);
	} else {
	germline.data <- NULL;
	}

if (length(somatic.files) > 0) {
	somatic.data <- read.delim(somatic.files[1]);
	} else {
	somatic.data <- NULL;
	}

# find and read in CNA data
cna.files <- list.files(path = '../CNAs', pattern = arguments$project, full.names = TRUE);
mops.files <- rev(sort(cna.files[grepl('mops__gene_data.tsv', cna.files)]));

if (length(mops.files) > 0) {
	cna.data <- read.delim(mops.files[1]);
	} else {
	cna.data <- NULL;
	}

# find and read in LOH data
loh.files <- list.files(path = '../SUMMARY', pattern = arguments$project, full.names = TRUE);
loh.files <- rev(sort(loh.files[grepl('LOH_data.tsv', loh.files)]));

if (length(loh.files) > 0) {
	loh.data <- read.delim(loh.files[1], check.names = FALSE);
	} else { 
	loh.data <- NULL;
	}

# find and read in Methylation data
methyl.files <- list.files(path = '../Methylation', pattern = arguments$project, full.names = TRUE);
gene.files <- rev(sort(methyl.files[grepl('gene_data.tsv', methyl.files)]));

if (length(gene.files) > 0) {
	methylation.data <- read.delim(gene.files[1]);
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
	keep.fields <- c('Sample','Hugo_Symbol','Region','Variant_Classification','HGVSc','HGVSp',
		'ClinVar','LOH','CN','t_vaf');

	somatic.data <- somatic.data[!grepl('Flank|UTR|Intron|Silent', somatic.data$Variant_Classification),];
	}

if (!is.null(germline.data)) {
	keep.fields <- c('Sample','Hugo_Symbol','Region','Variant_Classification','HGVSc','HGVSp',
		'ClinVar','LOH','CN','t_vaf');

	germline.data <- germline.data[grepl('pathogenic', germline.data$ClinVar),keep.fields];
	}

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
	tmp[,c('LOH','germline','somatic','methylation','copynumber')] <- NA;

	df1 <- unique(germline.data[which(germline.data$Hugo_Symbol == gene),]$Sample);
	df2 <- unique(somatic.data[which(somatic.data$Hugo_Symbol == gene),]$Sample);
	df3 <- unique(methylation.data[which(methylation.data[,gene] > 0.1),]$Sample);
	df4 <- names(which(apply(loh.data[which(loh.data$Gene == gene),6:ncol(loh.data)],2,
		function(i) { any(na.omit(i) == 'LOH') } )));

	df5.up <- cna.data[which(cna.data[,gene] > 2.8),]$Sample;
	df5.dn <- cna.data[which(cna.data[,gene] < 1.7),]$Sample;

	if (length(df1) > 0) {
		tmp[which(tmp$Sample %in% df1),]$germline <- 'germline';
		}
	if (length(df2) > 0) {
		tmp[which(tmp$Sample %in% df2),]$somatic <- 'somatic';
		}
	if (length(df3) > 0) {
		tmp[which(tmp$Sample %in% df3),]$methylation <- 'methylation';
		}
	if (length(df4) > 0) {
		tmp[which(tmp$Sample %in% df4),]$LOH <- 'LOH';
		}
	if (length(df5.up) > 0) {
		tmp[which(tmp$Sample %in% df5.up),]$copynumber <- 'cn.gain';
		}
	if (length(df5.dn) > 0) {
		tmp[which(tmp$Sample %in% df5.dn),]$copynumber <- 'cn.loss';
		}

	tmp[,gene] <- apply(tmp[,c('germline','somatic','methylation','copynumber','LOH')],1,
		function(i) { paste(na.omit(i), collapse = '/') });

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
