#!/usr/bin/env perl
### pughlab_multimmr_pipeline.pl ###################################################################
use AutoLoader 'AUTOLOAD';
use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use YAML qw(LoadFile);
use List::Util qw(any none);
use Data::Dumper;

my $cwd = dirname(__FILE__);
require "$cwd/scripts/utilities.pl";

####################################################################################################
# version       author	  	comment
# 1.0		sprokopec       script to run PughLab DNASeq + EMSeq pipeline

### USAGE ##########################################################################################
# pughlab_dnaseq_pipeline.pl -c tool_config.yaml -d data.yaml
#
# where:
#	- dnaseq_tool_config.yaml contains tool versions and parameters, output directory, reference
#	information, etc. for the DNA-seq pipeline
#	- dnaseq_data_config.yaml contains sample information (YAML file containing paths to FASTQ 
#	files) for the DNA-seq pipeline
#	- emseq_tool_config.yaml contains tool versions and parameters, output directory, reference
#	information, etc. for the EM-seq pipeline
#	- emseq_data_config.yaml contains sample information (YAML file containing paths to FASTQ 
#	files) for the EM-seq pipeline

### DEFINE SUBROUTINES #############################################################################
# run the DNA-seq pipeline 
sub run_dnaseq_pipeline {
	my %args = (
		tool	=> undef,
		data	=> undef,
		cluster	=> 'slurm',
		cleanup => 1,
		@_
		);

	my $command = join(' ',
		"perl $cwd/pughlab_dnaseq_pipeline.pl",
		"--tool", $args{tool},
		"--data", $args{data},
		"--preprocessing --qc --variant_calling",
		"-c", $args{cluster}
		);

	if ($args{cleanup}) {
		$command .= " --remove";
		}

	return($command);
	}

# run the EM-seq pipeline 
sub run_emseq_pipeline {
	my %args = (
		tool	=> undef,
		data	=> undef,
		cluster	=> 'slurm',
		cleanup => 1,
		@_
		);

	my $command = join(' ',
		"perl $cwd/pughlab_emseq_pipeline.pl",
		"--tool", $args{tool},
		"--data", $args{data},
		"--fastq_prep --alignment --qc --analysis",
		"-c", $args{cluster}
		);

	if ($args{cleanup}) {
		$command .= " --remove";
		}

	return($command);
	}

# run the report-generation pipeline 
sub run_report_pipeline {
	my %args = (
		wgs_tool	=> undef,
		wgs_data	=> undef,
		dna_tool	=> undef,
		dna_data	=> undef,
		em_tool		=> undef,
		em_data		=> undef,
		report_dir	=> undef,
		cluster		=> 'slurm',
		setup		=> 0,
		@_
		);

	my $command = join(' ',
		"perl $cwd/scripts/report/make_multimmr_reports.pl",
		"--wgs_tool", $args{wgs_tool},
		"--wgs_data", $args{wgs_data},
		"--dna_tool", $args{dna_tool},
		"--dna_data", $args{dna_data},
		"--em_tool", $args{em_tool},
		"--em_data", $args{em_data},
		"--output_directory", $args{report_dir},
		"-c", $args{cluster}
		);

	if ($args{setup}) {
		$command .= " --dry-run";
		}

	return($command);
	}

### MAIN ###########################################################################################
sub main {
	my %args = (
		wgs_tool_config		=> undef,
		wgs_data_config		=> undef,
		dna_tool_config		=> undef,
		dna_data_config		=> undef,
		em_tool_config		=> undef,
		em_data_config		=> undef,
		hpc_driver		=> undef,
		del_intermediates	=> undef,
		dry_run			=> undef,
		@_
		);

	### PREAMBLE ######################################################################################
	my $date = strftime "%F", localtime;
	my $timestamp = strftime "%F_%H-%M-%S", localtime;

	# load tool config
	my $tool_data = LoadFile($args{dna_tool_config});

	# organize output and log directories
	my $output_directory = dirname($tool_data->{output_dir});
	my $log_directory = join('/', $output_directory, 'logs', 'run_MultiMMR_pipelines_' . $timestamp);

	unless(-e $log_directory) { make_path($log_directory); }

	# start logging
	my $log_file = join('/', $log_directory, 'run_MultiMMR_pipeline.log');
	open (my $log, '>', $log_file) or die "Could not open $log_file for writing.";

	print $log "---\n";
	print $log "Running MultiMMR pipeline.\n";

	if ( (defined($args{wgs_tool_config})) & (defined($args{wgs_data_config})) ) {
		print $log "\n  Config used for sWGS: $args{wgs_tool_config}";
		print $log "\n  Sample config used for sWGS: $args{wgs_data_config}\n";
		}

	if ( (defined($args{dna_tool_config})) & (defined($args{dna_data_config})) ) {
		print $log "\n  Config used for DNA-Seq: $args{dna_tool_config}";
		print $log "\n  Sample config used for DNA-Seq: $args{dna_data_config}\n";
		}

	if ( (defined($args{em_tool_config})) & (defined($args{em_data_config})) ) {
		print $log "\n  Config used for EM-Seq: $args{em_tool_config}";
		print $log "\n  Sample config used for EM-Seq: $args{em_data_config}";
		}

	print $log "\n---\n\n";

	# set tools and versions
	my $perl	= 'perl/' . $tool_data->{perl_version};
	my $samtools	= 'samtools/' . $tool_data->{samtools_version};

        # indicate maximum time limit for parent jobs to wait
	my $max_time = '5-00:00:00';
	my $max_mem = '2G';

	### RUN ###########################################################################################
	# create objects to hold run ids
	my ($wgs_run_id, $dna_run_id, $em_run_id, $final_run_id);

	# create other objects
	my ($run_script);

	# is sWGS data provided?
	if (defined($args{wgs_tool_config}) & defined($args{wgs_data_config})) {

		# write command to run the DNA-Seq pipeline
		my $wgs_command = run_dnaseq_pipeline(
			tool	=> $args{wgs_tool_config},
			data	=> $args{wgs_data_config},
			cluster	=> $args{hpc_driver},
			cleanup => $args{del_intermediates}
			);

		# record command (in log directory) and then run job
		print $log "\n >> Submitting job for sWGS pipeline...\n";
		print $log "  COMMAND: $wgs_command\n\n";

		$run_script = write_script(
			log_dir	=> $log_directory,
			name	=> 'run_sWGS_pipeline',
			cmd	=> $wgs_command,
			modules	=> [$perl, $samtools],
			max_time	=> $max_time,
			mem		=> $max_mem,
			hpc_driver	=> $args{hpc_driver}
			);

		if ($args{dry_run}) {

			$wgs_command .= " --dry-run";
			`$wgs_command`;
			$wgs_run_id = 'run_sWGS_pipeline';

			} else {
			$wgs_run_id = submit_job(
				jobname		=> 'run_sWGS_pipeline',
				shell_command	=> $run_script,
				hpc_driver	=> $args{hpc_driver},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			print $log ">>> sWGS pipeline job id: $wgs_run_id\n\n";
			}
		}

	# is DNA-Seq (panel) data provided?
	if (defined($args{dna_tool_config}) & defined($args{dna_data_config})) {

		# write command to run the DNA-Seq pipeline
		my $dna_command = run_dnaseq_pipeline(
			tool	=> $args{dna_tool_config},
			data	=> $args{dna_data_config},
			cluster	=> $args{hpc_driver},
			cleanup => $args{del_intermediates}
			);

		# record command (in log directory) and then run job
		print $log "\n >> Submitting job for DNA-Seq pipeline...\n";
		print $log "  COMMAND: $dna_command\n\n";

		$run_script = write_script(
			log_dir	=> $log_directory,
			name	=> 'run_dnaseq_pipeline',
			cmd	=> $dna_command,
			modules	=> [$perl, $samtools],
			max_time	=> $max_time,
			mem		=> $max_mem,
			hpc_driver	=> $args{hpc_driver}
			);

		if ($args{dry_run}) {

			$dna_command .= " --dry-run";
			`$dna_command`;
			$dna_run_id = 'run_dnaseq_pipeline';

			} else {
			$dna_run_id = submit_job(
				jobname		=> 'run_dnaseq_pipeline',
				shell_command	=> $run_script,
				hpc_driver	=> $args{hpc_driver},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			print $log ">>> DNA-Seq pipeline job id: $dna_run_id\n\n";
			}
		}

	# is EM-Seq (panel) data provided?
	if (defined($args{em_tool_config}) & defined($args{em_data_config})) {

		# write command to run the EM-Seq pipeline
		my $em_command = run_emseq_pipeline(
			tool	=> $args{em_tool_config},
			data	=> $args{em_data_config},
			cluster	=> $args{hpc_driver},
			cleanup => $args{del_intermediates}
			);

		# record command (in log directory) and then run job
		print $log "\n >> Submitting job for EM-Seq pipeline...\n";
		print $log "  COMMAND: $em_command\n\n";

		$run_script = write_script(
			log_dir	=> $log_directory,
			name	=> 'run_emseq_pipeline',
			cmd	=> $em_command,
			modules	=> [$perl, $samtools],
			max_time	=> $max_time,
			mem		=> $max_mem,
			hpc_driver	=> $args{hpc_driver}
			);

		if ($args{dry_run}) {

			$em_command .= " --dry-run";
			`$em_command`;
			$em_run_id = 'run_emseq_pipeline';

			} else {
			$em_run_id = submit_job(
				jobname		=> 'run_emseq_pipeline',
				shell_command	=> $run_script,
				hpc_driver	=> $args{hpc_driver},
				dry_run		=> $args{dry_run},
				log_file	=> $log
				);

			print $log ">>> EM-Seq pipeline job id: $em_run_id\n\n";
			}
		}

	# write command to run the final report generation
	my $report_command = run_report_pipeline(
		wgs_tool	=> $args{wgs_tool_config},
		dna_tool	=> $args{dna_tool_config},
		em_tool		=> $args{em_tool_config},
		report_dir	=> $args{report_directory},
		cluster		=> $args{hpc_driver},
		setup		=> $args{dry_run}
		);

	# record command (in log directory) and then run job
	print $log "\n >> Submitting job for MultiMMR Report pipeline...\n";
	print $log "  COMMAND: $report_command\n\n";

	$run_script = write_script(
		log_dir	=> $log_directory,
		name	=> 'run_report_pipeline',
		cmd	=> $report_command,
		modules	=> [$perl, $samtools],
		dependencies	=> join(':', $wgs_run_id, $dna_run_id, $em_run_id),
		max_time	=> $max_time,
		mem		=> $max_mem,
		hpc_driver	=> $args{hpc_driver}
		);

	if ($args{dry_run}) {
		$final_run_id = 'run_report_pipeline';
		} else {
		$final_run_id = submit_job(
			jobname		=> 'run_report_pipeline',
			shell_command	=> $run_script,
			hpc_driver	=> $args{hpc_driver},
			dry_run		=> $args{dry_run},
			log_file	=> $log
			);

		print $log ">>> MultiMMR Report pipeline job id: $final_run_id\n\n";
		}


	# finish up
	print $log "\nProgramming terminated successfully.\n\n";
	close $log;

	}

### GETOPTS AND DEFAULT VALUES #####################################################################
# declare variables
my ($wgs_tool_config, $wgs_data_config);
my ($dna_tool_config, $dna_data_config, $em_tool_config, $em_data_config);
my ($dry_run, $help, $no_wait, $remove_intermediates, $report_directory);
my $hpc_driver = 'slurm';

# get command line arguments
GetOptions(
	'h|help'	=> \$help,
	'w|tool_wgs=s'	=> \$wgs_tool_config,
	't|tool_dna=s'	=> \$dna_tool_config,
	'f|tool_em=s'	=> \$em_tool_config,
	's|data_wgs=s'	=> \$wgs_data_config,
	'd|data_dna=s'	=> \$dna_data_config,
	'e|data_em=s'	=> \$em_data_config,
	'r|report_dir=s'	=> \$report_directory,
	'c|cluster=s'	=> \$hpc_driver,
	'remove'	=> \$remove_intermediates,
	'dry-run'	=> \$dry_run
	);

if ($help) {
	my $help_msg = join("\n",
		"Options:",
		"\t--help|-h\tPrint this help message",
		"\t--tool_wgs|-w\t<string> tool config for sWGS (yaml format)",
		"\t--data_wgs|-s\t<string> data config for sWGS (yaml format)",
		"\t--tool_dna|-t\t<string> tool config for DNA-Seq (yaml format)",
		"\t--data_dna|-d\t<string> data config for DNA-Seq (yaml format)",
		"\t--tool_em|-f\t<string> tool config for EM-Seq (yaml format)",
		"\t--data_em|-e\t<string> data config for EM-Seq (yaml format)",
		"\t--report_dir|-r\t<string> path to directory for patient reports",
		"\t--cluster|-c\t<string> cluster scheduler (default: slurm)",
		"\t--remove\t<boolean> should intermediates be removed? (default: false)",
		"\t--dry-run\t<boolean> should jobs be submitted? (default: false)"
		);

	print "$help_msg\n";
	exit;
	}

# do some quick error checks to confirm valid arguments	
if ( (!defined($wgs_tool_config)) | (!defined($wgs_data_config)) ) { 
	warn("One or more config files for sWGS are missing; will skip this portion.");
	}

if ( (!defined($dna_tool_config)) | (!defined($dna_data_config)) ) { 
	warn("One or more config files for DNA-Seq are missing; will skip this portion.");
	}

if ( (!defined($em_tool_config)) | (!defined($em_data_config)) ) { 
	warn("One or more config files for EM-Seq are missing; will skip this portion.");
	}

if (!defined($report_directory)) {
	die("No output directory defined for final patient reports; please provide --report_dir");
	}

main(
	wgs_tool_config		=> $wgs_tool_config,
	wgs_data_config		=> $wgs_data_config,
	dna_tool_config		=> $dna_tool_config,
	dna_data_config		=> $dna_data_config,
	em_tool_config		=> $em_tool_config,
	em_data_config		=> $em_data_config,
	report_directory	=> $report_directory,
	hpc_driver		=> $hpc_driver,
	del_intermediates	=> $remove_intermediates,
	dry_run			=> $dry_run
	);

