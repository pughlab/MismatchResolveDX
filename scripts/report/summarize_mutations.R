### summarize_mutations.R ##########################################################################
# Annotate each region with LOH status and then extract significant germline/somatic mutations.

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

# function to trim sample IDs
simplify.ids <- function(x) {
	match <- TRUE;
	if (length(x) == 1) {
		match <- FALSE;
		new.ids <- x;
		}
	index <- 1;
	while (match) {
		if (length(unique(sapply(x,function(i) { unlist(strsplit(i,''))[index] } ))) == 1) {
			index <- index + 1;
			} else {
			new.ids <- sapply(x,function(i) { substring(i,index,nchar(i)) } );
			match <- FALSE;
			}
		}
	return(new.ids);
	}

### PREPARE SESSION ################################################################################
# import command line arguments
library(argparse);

parser <- ArgumentParser();

parser$add_argument('-p', '--project', type = 'character', help = 'PROJECT name');
parser$add_argument('-o', '--output_dir', type = 'character', help = 'path to output directory');
parser$add_argument('-r', '--report_dir', type = 'character', help = 'path to report directory');
parser$add_argument('-y', '--sample_yaml', type = 'character', help = 'path to sample (BAM) yaml');
parser$add_argument('-m', '--mutations', type = 'character', help = 'path to ensemble mutation data');
parser$add_argument('-g', '--germline', type = 'character', help = 'path to germline mutation data');
parser$add_argument('-s', '--summary', type = 'character', help = 'path to summarized panel data');
parser$add_argument('-i', '--ichor', type = 'character', help = 'path to ichorCNA TF estimates');
parser$add_argument('--msi', type = 'character', help = 'path to MSI estimates');

arguments <- parser$parse_args();

###
#arguments$project <- 'UHNBB';
#arguments$output_dir <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/Patient_Reports/Mutations';
#arguments$report_dir <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/Patient_Reports/Reports';
#arguments$sample_yaml <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/UHNBB/DynaCare/DNASeq/GATK/gatk_bam_config_2025-08-25_15-25-35.yaml';
#arguments$mutations <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/Patient_Reports/data/DNASeq/ENSEMBLE/UHNBB_ensemble_mutation_data.tsv';
#arguments$germline <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/UHNBB/DynaCare/DNASeq/HaplotypeCaller/CPSR/2025-08-27_MultiMMR_mutations_for_cbioportal.tsv';
#arguments$summary <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/Patient_Reports/SUMMARY/UHNBB_summarized_panel_data.RData';
#arguments$ichor <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/UHNBB/DynaCare/sWGS/IchorCNA/2026-01-14_MultiMMR_ichorCNA_estimates.tsv';
#arguments$msi <- '/cluster/projects/pughlab/pughlab/projects/MultiMMR/UHNBB/DynaCare/DNASeq/panelMSI/2025-12-16_MultiMMR_panelMSI_estimates.tsv';
###


# load libraries
library(GenomicRanges);
library(plyr);
library(yaml);
library(xtable);

# what's the date?
date <- Sys.Date();

# get input files
output.dir <- arguments$output_dir;
reports.dir <- arguments$report_dir;

# create (if necessary) and move to output directory
if (!dir.exists(output.dir)) {
	dir.create(output.dir);
	}

if (!dir.exists(reports.dir)) {
	dir.create(reports.dir);
	}

genes.of.interest <- 'MLH1|MSH2|MSH3|MSH6|PMS2|APC|POLD1|POLE|MUTYH|TP53|BRAF';

### READ DATA ######################################################################################
# get data
if (is.null(arguments$mutations)) {
	stop('ERROR: No mutation provided, please provide path to variant MAF file.');
	} else {
	snp.data <- read.delim(arguments$mutations, skip = 1);
	}

if (!is.null(arguments$germline)) {
	germline.data <- read.delim(arguments$germline);
	} else {
	germline.data <- NULL;
	}

if (!is.null(arguments$summary)) {
	load(arguments$summary);
	} else {
	sample.data <- NULL;
	}

if (!is.null(arguments$ichor)) {
	ichor.data <- read.delim(arguments$ichor);
	} else {
	ichor.data <- NULL;
	}

if (!is.null(arguments$msi)) {
	msi.data <- read.delim(arguments$msi);
	} else {
	msi.data <- NULL;
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
tumour.smps <- sample.info[which(sample.info$Type == 'tumour'),]$Sample;

setwd(output.dir);


### FORMAT DATA ####################################################################################
# find LOH
loh.calls <- NULL;

if (!is.null(sample.data)) {

	loh.calls <- lapply(tumour.smps, function(smp) {

		if (smp %in% names(sample.data)) {

			smp.data <- sample.data[[smp]];
			
			# are all B-allele frequencies of germline heterozygous SNPs now homozygous?
			snp.loh <- merge(
				aggregate(BAF ~ Exon, smp.data, function(i) { all(i > 0.8) | all(i < 0.2) } ),
				aggregate(BAF ~ Exon, smp.data, function(i) { median(i, na.rm = TRUE) } ),
				by = 'Exon',
				all = TRUE
				);

			snp.loh$LOH.snp <- as.numeric(snp.loh$BAF.x) * sign(snp.loh$BAF.y - 0.5);

			# is the region copy-altered?
			cn.loh <- merge(
				aggregate(CN.Seg ~ Exon, smp.data, function(i) { all(i > 0.8) | all(i < -0.8) } ),
				aggregate(CN.Seg ~ Exon, smp.data, function(i) { median(i, na.rm = TRUE) } ),
				by = 'Exon',
				all = TRUE
				);

			cn.loh$LOH.cn <- as.numeric(cn.loh$CN.Seg.x) * sign(cn.loh$CN.Seg.y);

			# format LOH calls
			loh.data <- merge(
				snp.loh[,c('Exon','LOH.snp')],
				cn.loh[,c('Exon','LOH.cn')],
				by = 'Exon',
				all = TRUE
				);

			loh.data$LOH <- sapply(loh.data$LOH.snp, function(i) { 
				ifelse (is.na(i), 'n/a', ifelse( (abs(i) > 0), 'Yes', 'No'))
				});

			loh.data$CN <- sapply(loh.data$LOH.cn, function(i) {
				ifelse (is.na(i), 'n/a', 
					ifelse( (i > 0), 'CNgain', 
						ifelse( (i < 0), 'CNloss', 'CNneutral')))
				});
	
			# if TF is too low, don't trust LOH classifications
			if (smp %in% ichor.data[which(ichor.data$Tumour.Fraction < 0.1),]$ID) {
				loh.data$LOH <- 'n/a';
				}

			# if this is a tumour-only case (no matched normal), don't trust LOH classifications
			if (!any('normal' %in% sample.info[which(sample.info$Sample == smp),]$Type)) {
				loh.data$LOH <- 'n/a';
				}

			# extract region details
			loh.data$Chromosome <- sapply(loh.data$Exon, function(i) { 
				rev(unlist(strsplit(as.character(i),'\\.')))[3] } );
			loh.data$Start <- as.numeric(sapply(loh.data$Exon, function(i) { 
				rev(unlist(strsplit(as.character(i),'\\.')))[2] } ));
			loh.data$End <- as.numeric(sapply(loh.data$Exon, function(i) { 
				rev(unlist(strsplit(as.character(i),'\\.')))[1] } ));

			colnames(loh.data)[which(colnames(loh.data) == 'Exon')] <- 'Region';

			return(loh.data[,c('Region','Chromosome','Start','End','LOH','CN')]);
			} else { 
			return(NULL);
			}
		});

	names(loh.calls) <- tumour.smps;
	}

# combine somatic/germline data
snp.data$Mutation_Status <- 'somatic';
germline.data$Mutation_Status <- 'germline';

missing.cols <- setdiff(colnames(germline.data), colnames(snp.data));
if (length(missing.cols) > 0) { snp.data[,missing.cols] <- NA; }

missing.cols <- setdiff(colnames(snp.data), colnames(germline.data));
if (length(missing.cols) > 0) { germline.data[,missing.cols] <- NA; }

snp.data <- rbind(snp.data, germline.data[,colnames(snp.data)]);

# extract clinvar pathogenic hits
snp.data$ClinVar <- sapply(snp.data$CLIN_SIG, function(i) { 
	if (is.na(i)) { return(NA) 
		} else {
		tmp <- unlist(strsplit(i,','));
		x <- if (any(tmp == 'pathogenic')) { 'pathogenic'
			} else if (any(tmp == 'likely_pathogenic')) { 'likely_pathogenic'
			} else if (any(tmp %in% c('uncertain_significance','conflicting_interpretations_of_pathogenicity'))) { 'VUS'
			} else if (any(tmp == 'likely_benign')) { 'likely_benign'
			} else if (any(tmp == 'benign')) { 'benign'
			} else { NA }
		return(x);
		}
	});

# subset and format mutation data
keep.fields <- c('Tumor_Sample_Barcode','Hugo_Symbol','Chromosome','Start_Position','End_Position',
	'Variant_Classification','Variant_Type','Reference_Allele','Tumor_Seq_Allele2','Mutation_Status',
	'HGVSc','HGVSp_Short','dbSNP_RS','ClinVar','t_depth','t_alt_count','n_depth','n_alt_count');

mutation.data <- do.call(rbind, lapply(tumour.smps, function(smp) {

	if (smp %in% names(loh.calls)) {

		snp.gr <- makeGRangesFromDataFrame(
			snp.data[which(snp.data$Tumor_Sample_Barcode == smp),keep.fields],
			seqnames.field = 'Chromosome',
			start.field = 'Start_Position',
			end.field = 'End_Position',
			keep.extra.columns = TRUE
			);

		loh.gr <- makeGRangesFromDataFrame(
			loh.calls[[smp]],
			seqnames.field = 'Chromosome',
			start.field = 'Start',
			end.field = 'End',	
			keep.extra.columns = TRUE
			);

		tmp <- suppressWarnings(mergeByOverlaps(snp.gr, loh.gr));

		tmp2 <- data.frame(tmp$snp.gr, tmp$Region, tmp$LOH, tmp$CN);
		colnames(tmp2) <- gsub('tmp\\.', '', colnames(tmp2));
		colnames(tmp2)[which(colnames(tmp2) == 'seqnames')] <- 'Chromosome';
		colnames(tmp2)[which(colnames(tmp2) == 'start')] <- 'Start_Position';
		colnames(tmp2)[which(colnames(tmp2) == 'end')] <- 'End_Position';
	
		return(tmp2);

		} else {

		return(snp.data[,keep.fields]);

		}})
	);

colnames(mutation.data)[which(colnames(mutation.data) == 'HGVSp_Short')] <- 'HGVSp';

# calculate VAFs
mutation.data$t_vaf <- mutation.data$t_alt_count / mutation.data$t_depth;
mutation.data$n_vaf <- mutation.data$n_alt_count / mutation.data$n_depth;


# add patient annotation
annotated.data <- merge(
	sample.info[which(sample.info$Type == 'tumour'),1:2],
	mutation.data,
	by.x = 'Sample',
	by.y = 'Tumor_Sample_Barcode',
	all.y = TRUE
	);

write.table(
	annotated.data[which(annotated.data$Mutation_Status == 'germline'),],
	file = generate.filename(arguments$project, 'germline_mutations', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

write.table(
	annotated.data[which(annotated.data$Mutation_Status == 'somatic'),],
	file = generate.filename(arguments$project, 'somatic_mutations', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# filter out silent mutations
annotated.data <- annotated.data[which(!annotated.data$Variant_Classification %in% c(
	'Intron','Silent','IGR','RNA')),];

annotated.data <- annotated.data[which(annotated.data$t_vaf >= 0.05),];


# annotate MSI data
msi.data <- merge(
	sample.info[which(sample.info$Type == 'tumour'),1:2],
	msi.data,
	by = 'Sample',
	all = TRUE
	);

msi.data$Status <- 'MSS';
if (any(msi.data$Proportion > 0.1)) {
	msi.data[which(msi.data$Proportion > 0.1),]$Status <- 'MSI';
	}


### PER-SAMPLE SUMMARIES ###########################################################################
setwd(reports.dir);

germline.fields <- c('Sample','Hugo_Symbol','Chromosome','Start_Position','End_Position',
	'HGVSc','HGVSp','Variant_Classification','ClinVar','LOH','CN','n_vaf','t_vaf');
somatic.fields <- c('Sample','Hugo_Symbol','Chromosome','Start_Position','End_Position',
	'HGVSc','HGVSp','Variant_Classification','ClinVar','LOH','CN','t_vaf');

# initiate report object
tex.file <- 'mutation_results.tex';

for (patient in unique(sample.info$Patient)) {

	if (!dir.exists(patient)) {
		dir.create(patient);
		}

	setwd(patient);

	# get/print germline variants
	germ.mutations <- annotated.data[which(annotated.data$Patient == patient & 
		annotated.data$Mutation_Status == 'germline'),germline.fields];
	germ.mutations$t_vaf <- round(germ.mutations$t_vaf * 100, 3);
	germ.mutations$n_vaf <- round(germ.mutations$n_vaf * 100, 3);

	write.table(
		germ.mutations,
		file = generate.filename(patient, 'germline_mutations','tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);

	germ.mutations <- germ.mutations[grepl('pathogenic', germ.mutations$ClinVar),];

	caption <- ifelse (nrow(germ.mutations) > 0, 
		'Germline mutations identified using the HaplotypeCaller/CPSR method; mutations shown include known pathogenic variants with tumour VAF $\\geq$ 5\\%. If more then 10 variants remain after filtering, priority is given to variants in genes of interest.',
		'No germline mutations remain after filters (known pathogenic variants with VAF $\\geq$ 5\\%).'
		);

	if (nrow(germ.mutations) > 10) {
		roi.idx <- grep(genes.of.interest, germ.mutations$Hugo_Symbol);
		ns.idx <- grep('Nonsense', germ.mutations$Variant_Classification);
		ms.idx <- grep('Missense', germ.mutations$Variant_Classification);

		germ.mutations <- unique(rbind(
			germ.mutations[intersect(roi.idx,ns.idx),], germ.mutations[intersect(roi.idx,ms.idx),], germ.mutations[roi.idx,],
			germ.mutations));

		germ.mutations <- germ.mutations[1:10,];
		}

	# add data to tex file
	write("\\section{Germline Mutations}\n", file = tex.file);

	if (nrow(germ.mutations) > 0) {

		print(
			xtable(
				germ.mutations[,c(1,2,6,7,9,10,11)],
				caption = caption
				),
			file = tex.file,
			append = TRUE,
			include.rownames = FALSE,
			latex.environments = ""
			);

		} else {

		write(paste0(caption, "\n"), file = tex.file, append = TRUE);

		}

	# get/print somatic variants
	smp.mutations <- annotated.data[which(annotated.data$Patient == patient & 
		annotated.data$Mutation_Status == 'somatic'),somatic.fields];
	colnames(smp.mutations) <- gsub('t_vaf','VAF', colnames(smp.mutations));
	smp.mutations$VAF <- round(smp.mutations$VAF * 100, 3);

	write.table(
		smp.mutations,
		file = generate.filename(patient, 'somatic_mutations','tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);

	caption <- ifelse (nrow(smp.mutations) > 0, 
		'Somatic mutations identified using our ENSEMBLE method; mutations shown include non-silent variants with VAF $\\geq$ 5\\%. If more then 10 variants remain after filtering, priority is given to pathogenic variants and variants in genes of interest. LOH is only indicated if i) matched normal is present and ii) the variant overlaps a region with heterozygous germline SNPs.',
		'No somatic mutations remain after filters (non-silent variants with VAF $\\geq$ 5\\%).'
		);

	if (nrow(smp.mutations) > 10) {
		path.idx <- grep('pathogenic', smp.mutations$ClinVar);
		roi.idx <- grep(genes.of.interest, smp.mutations$Hugo_Symbol);
		ns.idx <- grep('Nonsense', smp.mutations$Variant_Classification);
		ms.idx <- grep('Missense', smp.mutations$Variant_Classification);

		smp.mutations <- unique(rbind(smp.mutations[path.idx,],
			smp.mutations[intersect(roi.idx,ns.idx),], smp.mutations[intersect(roi.idx,ms.idx),], smp.mutations[roi.idx,],
			smp.mutations));

		smp.mutations <- smp.mutations[1:10,];
		}

	# add data to tex file
	write("\\section{Somatic Mutations}\n", file = tex.file, append = TRUE);

	if (nrow(smp.mutations) > 0) {

		print(
			xtable(
				smp.mutations[,c(1,2,6,7,9,10,11)],
				caption = caption
				),
			file = tex.file,
			append = TRUE,
			include.rownames = FALSE,
			latex.environments = ""
			);

		} else {

		write(paste0(caption, "\n"), file = tex.file, append = TRUE);

		}

	# add MSI data to tex file
	write("\\subsection{MSI Status}\n", file = tex.file, append = TRUE);

	smp.msi <- msi.data[which(msi.data$Patient == patient),];
	caption <- if (nrow(smp.msi) > 0) { 'MSI status is based on proportion of target MS regions with reduced mononucleotide repeat count (compared to matched normal or panel of normals; proportion $\\geq$ 10\\%) is deemed unstable).'
		} else { 'MSI Status was not evaluated.'; }

	if (nrow(smp.msi) > 0) {

		print(
			xtable(
				smp.msi[,c('Sample','Status','Proportion')],
				caption = caption
				),
			file = tex.file,
			append = TRUE,
			include.rownames = FALSE,
			latex.environments = ""
			);

		} else {

		write(paste0(caption, "\n"), file = tex.file, append = TRUE);

		}

	write("\\pagebreak", file = tex.file, append = TRUE);

	setwd(reports.dir);
	}

### SAVE SESSION INFO ##############################################################################
setwd(output.dir);
save.session.profile(generate.filename(arguments$project, 'AnnotateLOH_SessionProfile','txt'));
