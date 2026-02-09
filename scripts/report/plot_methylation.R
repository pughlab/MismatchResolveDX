### plot_methylation.R #############################################################################
# Plot EM-Seq methylation estimates

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

# function to format data for output table
format.mean <- function(x) {
	paste0(round(median(x),2), ' (', round(quantile(x,0.25),2), ' - ', round(quantile(x,0.75),2), ')')
	}

### PREPARE SESSION ################################################################################
# import command line arguments
library(argparse);

parser <- ArgumentParser();

parser$add_argument('-p', '--project', type = 'character', help = 'project name');
parser$add_argument('-o', '--output_dir', type = 'character', help = 'path to output directory');
parser$add_argument('-r', '--report_dir', type = 'character', help = 'path to report directory');
parser$add_argument('-y', '--sample_yaml', type = 'character', help = 'path to sample (BAM) yaml');
parser$add_argument('-s', '--snp_data', type = 'character', help = 'path to ensemble mutations');
parser$add_argument('-m', '--methylation', type = 'character', help = 'path to promoter-level methylation estimates');


arguments <- parser$parse_args();

# load libraries
library(BoutrosLab.plotting.general);
library(yaml);
library(xtable);

# what's the date?
date <- Sys.Date();

# get input files
methyl.file <- arguments$methylation;
mutation.file <- arguments$snp_data;
output.dir <- arguments$output_dir;
reports.dir <- arguments$report_dir;

# create (if necessary) and move to output directory
if (!dir.exists(output.dir)) {
	dir.create(output.dir);
	}

if (!dir.exists(reports.dir)) {
	dir.create(reports.dir);
	}

genes.of.interest <- 'MLH1|MSH2|MSH3|MSH6|PMS2';

### MAIN ##########################################################################################
## read in data files
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
smp.names <- sample.info$Sample;


# read in methylation estimates
methylation.data <- read.delim(methyl.file, check.names = FALSE);
check.samples <- intersect(smp.names, colnames(methylation.data));

gene.data <- methylation.data[grepl(genes.of.interest, methylation.data$V4),];
rownames(gene.data) <- sapply(gene.data$V4, function(i) { unlist(strsplit(i,'_'))[1] } );

genes.of.interest <- unlist(strsplit(genes.of.interest,'\\|'));
gene.data <- gene.data[genes.of.interest,];

gene.data <- merge(
	sample.info,
	t(gene.data[,check.samples]),
	by.x = 'Sample',
	by.y = 'row.names',
	all.x = TRUE
	);

gene.data$Patient <- as.factor(gene.data$Patient);


# read in mutation data and extract BRAF V600E
if (!is.null(mutation.file)) {
	mutation.data <- read.delim(mutation.file, comment.char = '#');
	braf.samples <- unique(mutation.data[which(mutation.data$Hugo_Symbol == 'BRAF' & 
		mutation.data$HGVSp_Short == 'p.V600E'),]$Tumor_Sample_Barcode);
	rm(mutation.data);
	} else { 
	braf.samples <- c('none');
	}


# move to output directory
setwd(output.dir);

save(
	sample.info,
	methylation.data,
	gene.data,
	braf.samples,
	file = generate.filename(arguments$project, 'methylation_data','RData')
	);

# determine some parameters
axis.cex <- if (nrow(sample.info) <= 30) { 1
	} else if (nrow(sample.info) <= 80) { 0.8
	} else if (nrow(sample.info) <= 100) { 0.6
	} else if (nrow(sample.info) <= 140) { 0.5
	} else if (nrow(sample.info) <= 180) { 0.4
	} else { 0 };

# format for plotting
plot.data <- gene.data;
plot.data$Type <- factor(plot.data$Type, levels = c('normal','tumour'));
plot.data <- plot.data[order(plot.data$Type, plot.data$MLH1),];
plot.data$Sample <- factor(plot.data$Sample, levels = unique(plot.data$Sample));
plot.data$colour <- 'grey80';
plot.data[which(plot.data$Sample %in% braf.samples),]$colour <- 'red';

braf.key <- list(
	rect = list(col = c('red','grey80'), size = 1.5, border = 'black', height = 0.85),
	text = list(lab = c('BRAF V600E','BRAF wildtype'), cex = 1)
	);

plot.objects <- list();

for (gene in rev(genes.of.interest)) {

	plot.objects[[gene]] <- create.barplot(
		as.formula(paste(gene, '~ Sample | gene')),
		plot.data,
		col = plot.data$colour,
		abline.h = 0.1,
		abline.lty = 2
		);

	plot.objects[[gene]]$par.strip.text$lab <- gene;
	plot.objects[[gene]]$par.strip.text$fontface <- 'italic';
	}	

# combine them!
create.multiplot(
	plot.objects = plot.objects,
	xlab.label = NULL,
	ylab.label = 'Fraction of reads with methyation',
	ylab.padding = 3,
	ylab.cex = 1.5,
	xaxis.rot = 90,
	ylimits = c(0,1),
	yat = seq(0,1,0.2),
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0.5,0),
	xaxis.cex = axis.cex,
	yaxis.cex = 1,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	print.new.legend = TRUE,
	legend = list(
		inside = list(fun = draw.key, args = list(key = braf.key), x = 0.01, y = 0.96)
		),
	filename = generate.filename(arguments$project, 'methylation_estimates','png'),
	y.spacing = 2,
	height = 11,
	width = 11,
	resolution = 600
	);

### PER-SAMPLE SUMMARIES ###########################################################################
# find all sWGS files
all.files <- list.files(pattern = 'methylation_data.RData');

orig.file <- generate.filename(arguments$project, 'methylation_data','RData');
all.files <- setdiff(all.files, orig.file);

orig.sample.info <- sample.info;
orig.metric.data <- gene.data;

combined.metrics <- gene.data;

for (file in all.files) {
	load(file);
	combined.metrics <- unique(rbind(combined.metrics, gene.data));
	}

plot.data <- reshape(combined.metrics,
	direction = 'long',
	varying = list(4:ncol(combined.metrics)),
	v.names = 'Fraction',
	idvar = c('Patient','Sample','Type'),
	timevar = 'Gene',
	times = colnames(combined.metrics)[4:ncol(combined.metrics)]
	);

# get average scores across normal samples per gene
normal.data <- list();
for (gene in genes.of.interest) {
	normal.data[[gene]] <- mean(
		combined.metrics[which(combined.metrics$Type == 'normal'),gene], na.rm = TRUE);
	}

setwd(reports.dir);

# initiate report object
tex.file <- 'methylation_estimates_EMSeq.tex';

for (patient in unique(orig.sample.info$Patient)) {

	if (!dir.exists(patient)) {
		dir.create(patient);
		}

	these.smps <- orig.sample.info[which(orig.sample.info$Patient == patient),]$Sample;
	smp.data <- combined.metrics[which(combined.metrics$Sample %in% these.smps),];
	rownames(smp.data) <- smp.data$Sample;

	if (!('tumour' %in% smp.data$Type)) { next; }

	setwd(patient);

	# start with historical data
	to.write <- data.frame(rbind(
		apply(
			combined.metrics[which(combined.metrics$Type == 'normal'),genes.of.interest]*100,
			2,format.mean),
		apply(
			combined.metrics[which(combined.metrics$Type == 'tumour'),genes.of.interest]*100,
			2, format.mean)
		));

	row.names(to.write) <- c('all normals','all tumours');

	to.write <- rbind(to.write, round(smp.data[,genes.of.interest]*100,2));

	caption1 <- "Methylation level estimates (percent of reads with a converted C) by EMSeq:MethylDackel across all previously processed normal or tumour samples (reported as median + IQR) and for individual patient samples.";
	caption2 <- "Methylation level estimates (fraction of reads with a converted C) for each sample.";

	# add data to tex file
	write("\\section{Methylation}\n", file = tex.file);
	write(paste0("MMR genes tested: ", paste(genes.of.interest, collapse = ', '), "\n"),
		file = tex.file, append = TRUE);

	print(
		xtable(
			t(to.write),
			caption = caption1
			),
		file = tex.file,
		append = TRUE,
		include.rownames = TRUE,
		latex.environments = ""
		);

	setwd(reports.dir);
	}

### SAVE SESSION INFO ##############################################################################
setwd(output.dir);
save.session.profile(generate.filename(arguments$project, 'PlotMethylation_SessionProfile','txt'));
