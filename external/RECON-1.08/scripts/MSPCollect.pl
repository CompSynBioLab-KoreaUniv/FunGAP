#! /usr/bin/perl

if (!@ARGV) {
    die "usage: MSPCollect BLAST_output_file\n";
}

open (BLST, "$ARGV[0]") || die "usage: MSPCollect BLAST_output_file\nCan not open the BLAST_output_file $ARGV[0]\n";

$score_cutoff = 0;
$iden_cutoff = 0;

$score = -1;


while (<BLST>) {
    if (/^Query=\s+(\S+)/) {
	$query = $1;
    }
    elsif (/^>(\S+)/) {
	$sbjct_current = $1;
    }
    elsif (/^\s+Score =\s*(\S+)/) {
	if (&filter) {
	    &report;
	}
	$head_q_boolean = 1;
	$head_s_boolean = 1;
	$sbjct = $sbjct_current;
	$score_string = $1;
	if ($score_string =~ /(\S+)(e|E)\+(\d+)/) {
	    $score = $1 * exp($3*log(10));
	}
	else {$score = $score_string;}
    }
    elsif (/^\s+Identities =\s+(\d+)\/(\d+)/) {
	$iden = $1/$2*100;
    }
    elsif (/^Query:\s+(\d+)\s+(\S+)\s+(\d+)/) {
	if ($head_q_boolean) {
	    $head_query = $1;
	    $head_q_boolean = 0;
	}
	$tail_query = $3;
    }
    elsif (/^Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)/) {
	if ($head_s_boolean) {
	    $head_sbjct = $1;
	    $head_s_boolean = 0;
	}
	$tail_sbjct = $3;
    }
    elsif (/^Parameters:/) {
	if (&filter) {&report;}
	$score = -1;
    }
}


sub filter {
    if ($score > 0 && $score > $score_cutoff && $iden > $iden_cutoff) {return 1;}
    return 0;
}


sub report {
    printf("%06d %03d %05d %05d %s %05d %05d %s \n", $score, $iden, $head_query, $tail_query, $query, $head_sbjct, $tail_sbjct, $sbjct);
}
