#!/usr/bin/perl
# --------------------------------------------
# Alex Lomsadze
# Georgia Institute of Technology
# 2015
# 
# This program takes as input file with CDS coordinates in sorted GTF format
# (as in genemark.gtf) and calculates location of complete introns 
# from coding exons coordinates.
# 
# Input must be sorted!
# --------------------------------------------

use strict;
use warnings;
use Getopt::Long;

my $Version = "1.1";

# --------------------------------------------
my $in = '';
my $out = '';
my $v = '';
# --------------------------------------------
if ( $#ARGV == -1 ) { print PrintUsage($Version); exit 1; }

my $result = GetOptions
(
  'in=s'    => \$in,
  'out=s'   => \$out,
  'verbose' => \$v,
);

if ( $#ARGV != -1 ) { die "error, unexpected value found on cmd: @ARGV\n"; }
if ( !$result ) { die "\n"; }
if ( !$in or !$out ) { die "error in cmd: out or in is missing\n"; }
if ( $in eq $out ) { die "error in cmd: out equals in\n"; }

# --------------------------------------------
# global tmp variables
my $seq_id;
my $name;
my $type;
my $left;
my $right;
my $score;
my $strand;
my $phase;
my $attr;

my $gene_id;
my $transcript_id;
# --------------------------------------------

# parse file ones to get
# $h{ $transcript_id }{'n'} = number of CDS for each transcript ID
# parse file second time
# $h{ $transcript_id }{'c'} = keep here current CDS order (first, second, etc ...)

my %h;

ParseCodingGeneIDs( $in, \%h );

open( IN, $in ) || die "Can't open $in: $!\n";
open( OUT, ">", $out ) || die "Can't open $out: $!\n";

my $id = '';
my $last_pos = 0;

my $count_single_cds = 0;

while ( my $record = <IN> )
{
	next if ( $record =~ /^\s*$/ );
	next if ( $record =~ /^#/ );
	next if ( $record !~ /\t/ );
	ParseGFFLine( $record );
	next if ($type ne 'CDS' );
	ParseGTFAttr( $attr );
	
	# continue only for multy CDS genes
	#
	if ( $h{$transcript_id}{'n'} == 1 )
	{
		$count_single_cds += 1;
		next;
	}

	# initialize on first (initial or terminal) CDS
	# save ID and right boundary
	#
	if ( $h{$transcript_id}{'c'} == 0 )
	{
		if ( $id ) { print "error 1, unexpected: $id\n$record"; exit 1;}
		if ( $last_pos ) { print "error 2, unexpected: $last_pos\n$record"; exit 1;}
		
		$id = $transcript_id;
		$last_pos = $right;

		$h{$transcript_id}{'c'} = 1;
		next;
	}

	# check that CDS belongs to the same transcript
	#
	if ( $transcript_id ne $id ) { print "error 3, unexpected unsorted format: $record"; exit 1;}
	
	# intron is between saved (right) position of previous CDS and left boundary of current CDS
	#
	print OUT PrintIntron( $last_pos + 1, $left - 1 );

	$last_pos = $right;
	$h{$transcript_id}{'c'} += 1;

	# terminate on last (initial or terminal) CDS
	#
	if ( $h{$transcript_id}{'c'} == $h{$transcript_id}{'n'} )
	{
		$id = '';
		$last_pos = 0;
		next;
	}
}

close OUT;
close IN;

print "$count_single_cds  number of single CDS transcripts\n" if $v;

exit;

# --------------------------------------------
sub PrintIntron
{
	my $L = shift;
	my $R = shift;

	my $text = $seq_id  ."\t";
	$text .= ( $name ."\t" );
	$text .= "intron\t";
	$text .= ( $L ."\t" );
	$text .= ( $R ."\t" );
	$text .= ( 1 ."\t" );
	$text .= ( $strand ."\t" );
	$text .= ( $phase ."\t" );
	$text .= (  "gene_id \"". $gene_id ."\"; transcript_id \"". $transcript_id ."\";"  );

	return $text ."\n";
}
# --------------------------------------------
sub ParseCodingGeneIDs
{
	my $name = shift;
	my $ref = shift;

	my $cds_count = 0;
	my $transcript_count = 0;

	open( my $IN, $name ) || die "Can't open ifile $name: $!\n";
	while ( my $line = <$IN> )
	{
		next if ( $line =~ /^\s*$/ );
		next if ( $line =~ /^#/ );
		next if ( $line !~ /\t/ );
		ParseGFFLine( $line );
		next if ($type ne 'CDS' );
		ParseGTFAttr( $attr );
		
		if ( exists $ref->{$transcript_id}{'n'} )
		{
			$ref->{$transcript_id}{'n'} += 1;
		}
		else
		{
			$ref->{$transcript_id}{'n'} = 1;
			$ref->{$transcript_id}{'c'} = 0;
			$transcript_count += 1;
		}
			
		++$cds_count;
	}
	close $IN;
	
	print "$cds_count  number of CDS\n" if $v;
	print "$transcript_count  number of transcripts\n" if $v;
 }
# --------------------------------------------
sub ParseGFFLine
{
	my $line = shift;
	
	if ( $line =~/^(\S+)\t(.+?)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(.+)$/ )
	{
		$seq_id = $1;
		$name = $2;
		$type = $3;
		$left = $4;
		$right = $5;
		$score = $6;
		$strand = $7;
		$phase = $8;
		$attr = $9;
		
		$gene_id = '';
		$transcript_id = '';
	}
	else
	{
		print "error, GFF unexpected format: $line";
		exit 1;
	}
}
# --------------------------------------------
sub ParseGTFAttr
{
	my $line = shift;
	
	if ( $line =~ /gene_id \"(\S+?)\";/ )
	{
		$gene_id = $1;
	}
	
	if ( $line =~ /transcript_id \"(\S+?)\";/ )
	{
		$transcript_id = $1;
	}
	
	if ( !$gene_id or !$transcript_id )
	{
		print "error, unexpected format in attribute field: $line";
		exit 1;
	}
}
# --------------------------------------------
sub PrintUsage
{
	my $label = shift;
	
	my $txt = "Usage: $0  --in [file name]  --out [file name] 

  This program takes as input file with CDS coordinates in sorted GTF format (as in genemark.gtf) 
  and calculates location of complete introns from coding exons coordinates

  --v
version $label
";
	return $txt;
}
# --------------------------------------------
