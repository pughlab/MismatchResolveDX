### create_fastq_yaml.pl ##########################################################################
use AutoLoader 'AUTOLOAD';
use strict;
use warnings;
use Carp;
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use List::Util qw(any);
use YAML qw(DumpFile);
use Data::Dumper;

####################################################################################################
# version       author          comment
# 1.0		sprokopec	create config file in yaml format

### USAGE ##########################################################################################
# create_fastq_yaml.pl -i INPUT -o /path/to/OUTPUT_CONFIG.yaml
#
# where INPUT is a tab separated file containing sample info for use as YAML entries:
#
#	Patient | Sample | Type | Library | Lane | R1 | R2
#	OTB-001 | OTB-001-T | tumour | DYNA-OTB-001-T-DNA-C_S241 | L001 | /path/to/DYNA-OTB-001-T-DNA-C_S241_L001_R1.fastq.gz | /path/to/DYNA-OTB-001-T-DNA-C_S241_L001_R2.fastq.gz
#	OTB-001 | OTB-001-T | tumour | DYNA-OTB-001-T-DNA-C_S241 | L002 | /path/to/DYNA-OTB-001-T-DNA-C_S241_L002_R1.fastq.gz | /path/to/DYNA-OTB-001-T-DNA-C_S241_L002_R2.fastq.gz
#	OTB-001 | OTB-001-N | normal | DYNA-OTB-001-N-DNA-C_S15 | L001 | /path/to/DYNA-OTB-001-N-DNA-C_S15_L001_R1.fastq.gz | /path/to/DYNA-OTB-001-N-DNA-C_S15_L001_R2.fastq.gz

### GETOPTS PLUS ERROR CHECKING AND DEFAULT VALUES #################################################
# declare variables
my $sample_info;
my $output_config;

# read in command line arguments
GetOptions(
	'i|input=s'		=> \$sample_info,
	'o|output_config=s'	=> \$output_config
	);

# check for and remove trailing / from data dir
$data_directory =~ s/\/$//;

### HANDLING FILES #################################################################################
# open sample file
open(my $INPUT, '<', $sample_info) or die "Cannot open '$sample_info' file\n";
my $header = <$INPUT>;

my $smp_data = {};

# process each sample in SAMPLES
while (my $line = <$INPUT>) {

	chomp($line);
	
	# split line info
	my @file_info = split(/\t/, $line);
	my $subject = $file_info[0];
	my $sample = $file_info[1];
	my $type = $file_info[2];
	my $lib = $file_info[3];
	my $lane = $file_info[4];
	my $r1 = $file_info[5];
	my $r2 = $file_info[6];

	$smp_data->{$subject}->{$sample}->{type} = $type;
	$smp_data->{$subject}->{$sample}->{libraries}->{$lib}->{runlanes}->{$lane}->{fastq}->{R1} = $r1;
	$smp_data->{$subject}->{$sample}->{libraries}->{$lib}->{runlanes}->{$lane}->{fastq}->{R2} = $r2;

	}

close $INPUT;

local $YAML::Indent = 4;
DumpFile($output_config, $smp_data);


exit;
