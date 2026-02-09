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
parser$add_argument('-s', '--sample_yaml', type = 'character', help = 'path to sample (BAM) yaml');
parser$add_argument('-m', '--msi', type = 'character', help = 'path to MSI estimates');

arguments <- parser$parse_args();

# load libraries
library(BoutrosLab.plotting.general);
library(yaml);
library(xtable);

# what's the date?
date <- Sys.Date();

# get input files
msi.file <- arguments$msi;
output.dir <- arguments$output_dir;
reports.dir <- arguments$report_dir;

# create (if necessary) and move to output directory
if (!dir.exists(output.dir)) {
	dir.create(output.dir);
	}

if (!dir.exists(reports.dir)) {
	dir.create(reports.dir);
	}

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


# read in MSI estimates
msi.data <- merge(
	sample.info[which(sample.info$Type != 'normal'),],
	read.delim(msi.file),
	by = 'Sample',
	all = TRUE
	);

# move to output directory
setwd(output.dir);

# determine some parameters
axis.cex <- if (nrow(msi.data) <= 30) { 1
	} else if (nrow(msi.data) <= 80) { 0.8
	} else if (nrow(msi.data) <= 100) { 0.6
	} else if (nrow(msi.data) <= 140) { 0.5
	} else if (nrow(msi.data) <= 180) { 0.4
	} else { 0 };

#type.colours <- c('navajowhite2','#F8B4E3','lightseagreen');
#names(type.colours) <- c('normal','plasma','tumour');

type.colours <- c('navajowhite2','lightseagreen');
names(type.colours) <- c('normal','tumour');

my.legend <- legend.grob(
	legend = list(
		legend = list(
			colours = type.colours[1],
			labels = names(type.colours)[1]
			),
		legend = list(
			colours = type.colours[2],
			labels = names(type.colours)[2]
#			),
#		legend = list(
#			colours = type.colours[3],
#			labels = names(type.colours)[3]
			)
		),
	layout = c(2,1),
	size = 1.5
	);


tf.key <- list(
	rect = list(col = c('red','grey80'), size = 1.5, border = 'black', height = 0.85),
	text = list(lab = c(expression(''>='5%'),expression(''<'5%')), cex = 1),
	title = 'Tumour Fraction',
	cex.title = 1.1,
	x = 0.8,
	y = 0.95
	);

plot.data <- msi.data[order(msi.data$Proportion, decreasing = TRUE),];
plot.data$colour <- 'grey80';
#plot.data[which(plot.data$Sample %in% braf),]$colour <- 'red';
plot.data$Sample <- factor(plot.data$Sample, levels = unique(plot.data$Sample));

msi.plot <- create.barplot(
	Proportion ~ Sample,
	plot.data,
	col = plot.data$colour,
	ylimits = c(0,1),
	yat = seq(0,1,0.2),
	xlab.label = NULL,
	ylab.label = 'Proportion of MS sites\nthat are unstable',
	ylab.cex = 1.5,
	ylab.axis.padding = 2,
	xaxis.rot = 90,
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0.5,0),
	xaxis.cex = axis.cex,
	yaxis.cex = axis.cex,
#	abline.h = 0.1,
#	abline.lty = 2,
	key = tf.key,
	style = 'Nature'
	);

write.plot(
	msi.plot,		
	height = 4,
	width = 11,
	resolution = 200,
	filename = generate.filename(arguments$project, 'MSI_estimates', 'png')
	);

save(
	sample.info,
	methylation.data,
	gene.data,
	file = generate.filename(arguments$project, 'methylation_data','RData')
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

	# calculate difference in methylation between T and N
	avg.normal <- if (nrow(smp.data[which(smp.data$Type == 'normal'),]) > 1) {
		apply(smp.data[which(smp.data$Type == 'normal'),genes.of.interest],2,median, na.rm = TRUE);
		} else {
		as.numeric(smp.data[which(smp.data$Type == 'normal'),genes.of.interest]);
		}

	diff.data <- if (('normal' %in% smp.data$Type) & ('tumour' %in% smp.data$Type)) {
		smp.data[which(smp.data$Type == 'tumour'),genes.of.interest] - avg.normal;
		} else {
		smp.data[which(smp.data$Type == 'tumour'),genes.of.interest] -
		unlist(normal.data)[genes.of.interest];
		}

	# is this difference sufficiently large? [use a change of 10% as a threshold
	diff.data[which(abs(diff.data) < 0.1)] <- 0
	diff.data <- sign(diff.data);

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

	row.names(to.write) <- c('Average normal','Average tumour');

	to.write <- rbind(to.write, round(smp.data[,genes.of.interest]*100,2));

	for (smp in rownames(diff.data)) {
		to.write[paste0(smp,'_status'),] <- as.character(factor(
			as.numeric(diff.data[smp,]),
			levels = c(-1,0,1),
			labels = c('hypomethylation','no change','hypermethylation')
			));
		}

	caption1 <- "Methylation level estimates (percent of reads with a converted C) by EMSeq:MethylDackel across all previously processed normal or tumour samples (reported as median + IQR)."
	caption2 <- "Methylation level estimates (percent of reads with a converted C) by EMSeq:MethylDackel for this patient. If no matched normal is available, methylation status for each tumour is determined using historical average across all normals; hypomethylation = decrease $\\geq$ 10\\%; hypermethylation = increase $\\geq$ 10\\% relative to normal(s).";

	# plot per-sample methylation estimates
	smp.data$colour <- type.colours[match(smp.data$Type, names(type.colours))];
	smp.data$pch <- 21;
	if (nrow(smp.data[which(smp.data$Type == 'normal'),]) %in% 1:5) {
		smp.data[which(smp.data$Type == 'normal'),]$pch <- c(21,22,23,24,25)[1:nrow(smp.data[which(smp.data$Type == 'normal'),])];
		} 
	if (nrow(smp.data[which(smp.data$Type == 'tumour'),]) %in% 1:5) {
		smp.data[which(smp.data$Type == 'tumour'),]$pch <- c(21,22,23,24,25)[1:nrow(smp.data[which(smp.data$Type == 'tumour'),])];
		} 

	smp.key <- list(
		points = list(col = 'black',fill = smp.data$colour, cex = 1, pch = smp.data$pch),
		text = list(lab = smp.data$Sample, cex = 1)
		);

	plot.objects <- list();
	for (gene in genes.of.interest) {

		plot.objects[[gene]] <- make.violin.plot(
			data = plot.data,
			extra.data = smp.data,
			gene = gene
			);
		}

	# combine them!
	create.multiplot(
		plot.objects = plot.objects,
		plot.layout = c(3,2),
		xlab.label = NULL,
		ylab.label = expression('Percent Methylation'),
		ylab.padding = 3,
		ylab.cex = 1.5,
		xaxis.rot = 45,
		ylimits = c(0,100),
		yat = seq(0,100,20),
		xaxis.tck = c(0.5,0),
		yaxis.tck = c(0.5,0),
		xaxis.cex = 1,
		yaxis.cex = 1,
		xaxis.fontface = 'plain',
		yaxis.fontface = 'plain',
		print.new.legend = TRUE,
		legend = list(
			inside = list(fun = my.legend, x = 0.7, y = 0.95),
			inside = list(fun = draw.key, args = list(key = smp.key), x = 0.7, y = 0.85)
			),
		filename = generate.filename(arguments$project, 'methylation_estimates','png'),
		height = 5,
		width = 7,
		resolution = 200
		);

	caption3 <- "Violin plots depect methylation level estimates (percent of reads with a converted C) for all previously processed normal or tumour samples, while individual points show results for this patient.";

	# add data to tex file
	write("\\section{Methylation}\n", file = tex.file);
	write(paste0("MMR genes tested: ", paste(genes.of.interest, collapse = ', '), "\n"),
		file = tex.file, append = TRUE);

	print(
		xtable(
			t(to.write[1:2,]),
			caption = caption1
			),
		file = tex.file,
		append = TRUE,
		include.rownames = TRUE,
		latex.environments = ""
		);

	normal.smps <- smp.data[which(smp.data$Type == 'normal'),]$Sample;
	tumour.smps <- smp.data[which(smp.data$Type == 'tumour'),]$Sample;

	for (smp in tumour.smps) {

		print(
			xtable(
				t(to.write[c(normal.smps, smp, paste0(smp,'_status')),]),
				caption = caption2
				),
			file = tex.file,
			append = TRUE,
			include.rownames = TRUE,
			latex.environments = ""
			);
		}
	
	write("\\pagebreak\n", file = tex.file, append = TRUE);
	write("\\begin{figure}[h!]", file = tex.file, append = TRUE);
	write("\\begin{center}", file = tex.file, append = TRUE);
	write(paste0(
		"\\includegraphics[width=0.9\\textwidth]{",
		generate.filename(arguments$project, 'methylation_estimates','png'), '}'
		), file = tex.file, append = TRUE);
	write("\\end{center}", file = tex.file, append = TRUE);
	write(paste0(
		"\\caption{", caption2, "}"
		), file = tex.file, append = TRUE);
	write("\\end{figure}\n", file = tex.file, append = TRUE);

	setwd(reports.dir);
	}

### SAVE SESSION INFO ##############################################################################
setwd(output.dir);
save.session.profile(generate.filename(arguments$project, 'PlotMethylation_SessionProfile','txt'));
