### plot_cna_summary.R ############################################################################
# Plot CNA landscape (from ichorCNA or panelCN.mops).

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

	# write key variables to file
	cat("\n### VARIABLES #################################################################\n");
	cat(paste0('Input: ', arguments$cnas));
	cat(paste0('CNA type: ', arguments$type));

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

parser$add_argument('-p', '--project', type = 'character', help = 'project name');
parser$add_argument('-o', '--output_dir', type = 'character', help = 'path to output directory');
parser$add_argument('-r', '--report_dir', type = 'character', help = 'path to report directory');
parser$add_argument('-s', '--sample_yaml', type = 'character', help = 'path to sample (BAM) yaml');
parser$add_argument('-c', '--cnas', type = 'character', help = 'path to gene by sample matrix of cna calls');
parser$add_argument('-t', '--type', type = 'character', help = 'one of ichor or mops');
parser$add_argument('-a', '--annotation', type = 'character', help = 'path to cytoband locations (for use with ichorCNA only)');

arguments <- parser$parse_args();

# import libraries
library(BoutrosLab.plotting.general);
library(GenomicRanges);
library(yaml);

# indicate key genes
classic.genes <- unlist(strsplit('MLH1|MSH2|MSH6|PMS2','\\|'));
ancillary.genes <- unlist(strsplit('EPCAM|MLH3|POLD1|POLE|MUTYH|TP53|BRAF','\\|'));

genes.of.interest <- c(classic.genes, ancillary.genes);


### READ DATA ######################################################################################
# check data type
if (is.null(arguments$type)) {
	stop('ERROR: No data type provided; please indicate if data are from sWGS (ichorCNA) or targeted-panel (panelCN.mops).');
	} else if (!grepl('ichor|mops', tolower(arguments$type))) {
	stop('ERROR: Unrecognized data type; please indicate if data are from ichor or mops.');
	}

# get data
if (grepl('ichor', tolower(arguments$type))) {
	if (is.null(arguments$annotation)) {
		stop('ERROR: No annotation file provided; please provide path to cytoband file.');
		}
	anno.data <- read.delim(arguments$annotation);
	}

if (is.null(arguments$cnas)) {
	stop('ERROR: No CNA calls provided, please provide path to CNA matrix from IchorCNA or panelCN.mops.');
	} else {
	input.data <- read.delim(arguments$cnas);
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
sample.info$Type <- factor(sample.info$Type, levels = c('tumour','normal'));

# we only care about tumour samples, so only show those
sample.info <- sample.info[which(sample.info$Type == 'tumour'),];
all.samples <- sample.info$Sample;

# create (if necessary) and move to output directory
if (!dir.exists(arguments$output_dir)) {
	dir.create(arguments$output_dir);
	}

setwd(arguments$output_dir);

### FORMAT DATA ####################################################################################
# makes sure input.data is formatted
if (grepl('ichor', tolower(arguments$type))) {

	# get arm-level annotations
	anno.data$arm <- apply(anno.data[,c('chrom','cytoband')], 1, function(i) { 
		paste0(sub('chr','', i[1]), substr(i[2],0,1));
		});
	anno.gr <- makeGRangesFromDataFrame(
		unique(anno.data[,c('chrom','start','end','arm')]),
		seqnames.field = 'chrom',
		start.field = 'start',
		end.field = 'end',
		keep.extra.columns = TRUE
		);

	# initiate arm x sample matrix
	plot.data <- unique(anno.data[,c('chrom','arm')]);
	plot.data$chrom <- factor(plot.data$chrom, levels = paste0('chr',c(1:22,'X')));
	plot.data <- plot.data[order(plot.data$chrom, plot.data$arm),];
	rownames(plot.data) <- paste0('chr',plot.data$arm);

	# calculate weighted mean for each sample x arm
	for (smp in unique(input.data$ID)) {
		
		plot.data[,smp] <- NA;

		ichor.gr <- makeGRangesFromDataFrame(
			input.data[which(input.data$ID == smp),c('chrom','start','end','logR_Copy_Number')],
			seqnames.field = 'chrom',
			start.field = 'start',
			end.field = 'end',
			keep.extra.columns = TRUE
			);

		x <- mergeByOverlaps(ichor.gr, anno.gr);
		test.data <- unique(data.frame(x$anno.gr, x$logR_Copy_Number));

		for (arm in plot.data$arm) {

			tmp <- test.data[which(test.data$arm == arm),];

			plot.data[which(plot.data$arm == arm),smp] <- weighted.mean(
				x = tmp$x.logR_Copy_Number,
				w = tmp$width / sum(tmp$width),
				na.rm = TRUE
				) - 2;
			}
		}

	rm(tmp, test.data, x, anno.gr, ichor.gr);

	anno.data <- plot.data[,1:2];
	plot.data <- plot.data[,3:ncol(plot.data)];

	# format panelCN.mops data
	} else if (grepl('mops', tolower(arguments$type))) {

	# remove low quality calls
	input.data <- input.data[which(input.data$lowQual != 'lowQual'),];

	# format seqnames
	input.data$Chr <- factor(input.data$Chr, levels = paste0('chr',c(1:22)));

	# format CN values
	input.data$CN <- as.numeric(gsub('CN','', input.data$CN));
	input.data$LRR <- log2(input.data$RC.norm / input.data$medRC.norm);

	# indicate region of EPCAM-MSH2 deletion
	for (smp in unique(input.data$Sample)) {
		idx1 <- which(input.data$Sample == smp & input.data$Gene == 'EPCAM')[1];
		idx2 <- which(input.data$Sample == smp & input.data$Gene == 'MSH2')[1];

		input.data[idx1:idx2,]$Gene <- 'EPCAM';
		}

	# exclude 'other' sites (these are generally MSI sites)
	input.data <- input.data[!grepl('other', input.data$Exon),];

	# reshape to region x sample matix
	cna.data <- reshape(
		input.data[,c('Chr','Start','End','Gene','Exon','Sample','LRR')],
		direction = 'wide',
		timevar = 'Sample',
		idvar = c('Chr','Start','End','Gene','Exon')
		);

	colnames(cna.data) <- gsub('LRR.','',colnames(cna.data));
	rownames(cna.data) <- paste0('Region',1:nrow(cna.data));

	# trim out regions with mostly lowQual data (~5%)
	region.counts <- apply(cna.data[,6:ncol(cna.data)],1,function(i) { length(i[!is.na(i)])/length(i) } );
	to.keep <- which(region.counts > 0.05);
	cna.data <- cna.data[to.keep,];

	anno.data <- cna.data[order(cna.data$Chr, cna.data$Start),1:5];
	colnames(anno.data)[which(colnames(anno.data) == 'Chr')] <- 'chrom';
	plot.data <- cna.data[order(cna.data$Chr, cna.data$Start),6:ncol(cna.data)];

	# summarize gene-level CN status using weighted mean
	gene.cnas <- aggregate(
		CN ~ Sample + Gene,
		input.data,
		mean,
		na.rm = TRUE
		);

	gene.data <- reshape(
		gene.cnas[which(gene.cnas$Gene %in% genes.of.interest),],
		direction = 'wide',
		idvar = 'Sample',
		timevar = 'Gene'
		);
	colnames(gene.data) <- gsub('CN.','', colnames(gene.data));

	}

# check for missing samples
if (length(setdiff(all.samples,colnames(plot.data))) > 0) {
	missing.samples <- setdiff(all.samples,colnames(plot.data));
	plot.data[,missing.samples] <- NA;
	}

if (length(all.samples) == 1) {
	plot.data <- data.frame(cbind(plot.data, plot.data));
	colnames(plot.data) <- rep(all.samples,2);
	} else {
	plot.data <- plot.data[,all.samples];
	}

# remove bins with 0 data
#na.counts <- apply(cna.data,1,function(i) { length(i[!is.na(i)]); });
#cna.data <- cna.data[which(na.counts > 0),]

### PLOT DATA ######################################################################################
# grab some parameters
axis.cex <- if (ncol(plot.data) <= 30) { 1
	} else if (ncol(plot.data) <= 50) { 0.75
	} else if (ncol(plot.data) <= 80) { 0.5
	} else if (ncol(plot.data) <= 100) { 0
	} else { 0 };

chromosomes <- gsub('chr','',tolower(anno.data$chrom));
chr.breaks <- get.line.breaks(chromosomes);
covariate.colours <- force.colour.scheme(chromosomes, scheme = 'chromosome');

covariate.legends <- legend.grob(
	legends =  list(
		legend = list(
			colours = force.colour.scheme(c(1:22,'x'), scheme = 'chromosome'),
			labels = c(1:22,'X')
			)
		),
	size = 2
	);

# indicate output filename
outfile <- if (tolower(arguments$type) == 'ichor') {
	generate.filename(arguments$project, 'sWGS_ichorCNA_landscape','png');
	} else {
	generate.filename(arguments$project, 'panelCNmops_landscape','png');
	}

# create heatmap
create.heatmap(
	plot.data,
	cluster.dimensions = 'none',
	covariates.top = list(
		rect = list(
			col = 'transparent',
			fill = covariate.colours,
			lwd = 0
			)
		),
	covariates.top.grid.border = list(col = 'black', lwd = 1),
	covariates.top.grid.col = list(col = 'black', lwd = 1),
	covariates.top.col.lines = chr.breaks-0.5,
	inside.legend = list(fun = covariate.legends, x = 1.02, y = 1),
	yaxis.lab = if (length(all.samples) == 1) { all.samples } else {
		gsub('\\.','-', simplify.ids(colnames(plot.data))) },
	yat = if (length(all.samples) == 1) { 1.5 } else { TRUE },
	xaxis.lab = rep('',nrow(plot.data)),
	yaxis.cex = axis.cex,
	xaxis.tck = 0,
	yaxis.tck = if (axis.cex == 0) { 0 } else { 0.2 },
	yaxis.fontface = 'plain',
	axes.lwd = 1,
	col.colour = 'black',
	grid.col = TRUE,
	force.grid.col = TRUE,
	col.lines = chr.breaks,
	print.colour.key = TRUE,
	fill.colour = 'grey50',
	at = seq(-2,2,0.1),
	colour.scheme = c('blue','white','red'),
	colourkey.labels.at = seq(-2,2,1),
	colourkey.labels = seq(-2,2,1),
	colourkey.cex = 1,
	height = if (length(all.samples) > 12) { 7 } else { 5 },
	width = 12,
	resolution = 600,
	right.padding = 12,
	filename = outfile
	);

# write summary data to file
outfile2 <- if (tolower(arguments$type) == 'ichor') {
	generate.filename(arguments$project, 'sWGS_ichorCNA__arm-level','tsv');
	} else {
	generate.filename(arguments$project, 'panelCNmops__gene_data','tsv');
	}

if (tolower(arguments$type) == 'ichor') {
	to.write <- data.frame(Sample = colnames(plot.data), t(plot.data));
	} else {
	to.write <- gene.data[,c('Sample',genes.of.interest)];
	}

write.table(
	to.write,
	file = outfile2,
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### SAVE SESSION INFO ##############################################################################
if (tolower(arguments$type) == 'ichor') {
	save.session.profile(generate.filename('IchorCNA','SessionProfile','txt'));
	} else {
	save.session.profile(generate.filename('panelCNmops','SessionProfile','txt'));
	}
