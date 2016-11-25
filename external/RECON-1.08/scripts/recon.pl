#! /usr/bin/perl

$path = "";

if (!$path) {
    die "Please modifiy this script to specify the path to the binaries.\n";
}

if (@ARGV < 2) {
    die "usage: recon seq_name_list_file MSP_file integer\nSee 00README for details.\n";
}

open(SEQ, "$ARGV[0]") || die "usage: recon seq_name_list_file MSP_file integer\nCan not open the seq_name_list_file $ARGV[0].\nSee 00README for details.\n";
close SEQ;

open(MSP, "$ARGV[1]") || die "usage: recon seq_name_list_file MSP_file integer\nCan not open the MSP_file $ARGV[1].\nSee 00README for details.\n";
close MSP;


`rm -f core`;

`rm -rf summary`;
`mkdir summary`;



# prepare images

`rm -rf images`;
`mkdir images`;

`$path/imagespread $ARGV[0] $ARGV[1] $ARGV[2]`;

if ($?) {die "imagespread failed. Exit code $?\n";}

if (@ARGV < 3) {
    $sect = 1;
} else {
    $sect = $ARGV[2];
}

for ($i=1; $i<=$sect; $i++) {
   $spread = "images/spread" . $i;   
   `sort -k 3,3 -k 4n,4n -k 5nr,5nr $spread >> images/images_sorted`;
   if ($?) {die "sort failed for $spread.\n";}
}

`rm -f images/spread*`;



# initial definition of elements

`rm -rf ele_def_res`;
`mkdir ele_def_res`;

`$path/eledef $ARGV[0] $ARGV[1] single`;
if ($?) {die "eledef failed. Exit code $?\n";}



# re-defining elements

`rm -rf ele_redef_res`;
`mkdir ele_redef_res`;

`rm -f tmp tmp2`;
`ln -s ele_def_res tmp`;
`ln -s ele_redef_res tmp2`;

`$path/eleredef $ARGV[0]`;
if ($?) {die "eleredef failed. Exit code $?\n";}

`rm -f tmp tmp2`;



# re-defining edges

`rm -rf edge_redef_res`;
`mkdir edge_redef_res`;

`rm -f tmp tmp2`;
`ln -s ele_redef_res tmp`;
`ln -s edge_redef_res tmp2`;

`$path/edgeredef $ARGV[0]`;
if ($?) {die "edgeredef failed. Exit code $?\n";}



# famdef

`rm -f tmp tmp2`;
`ln -s edge_redef_res tmp`;

`$path/famdef $ARGV[0]`;
if ($?) {die "famdef failed. Exit code $?\n";}

`rm -f tmp`;

