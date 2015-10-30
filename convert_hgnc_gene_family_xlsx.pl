#!/usr/bin/env perl 
use Carp;
use LWP::Simple;
use File::Slurp;
use ConciseDescriptions;
#use Spreadsheet::WriteExcel;
use Excel::Writer::XLSX;
use strict;
use warnings;
#
# Authors: J. Done and R. Kishore, California Institute of Technology, 2014. 
# http://wiki.wormbase.org/index.php/Generation_of_automated_descriptions
#
# html path is defined by reading the an input file (e.g., docroot in RH, /var/www in Ubuntu, etc.)
my @args;
my $status;
#
my $species = "c_elegans";

my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
my $RELEASE = ConciseDescriptions::get_release();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
#
my $orthology_dir = $home . "release/$PRODUCTION_RELEASE/$species/orthology/input_files/";
my $hgnc_family_file  = $orthology_dir . "human_gene_families.txt";
my $xls_file = $hgnc_family_file;
   $xls_file =~s/\.txt/\.xlsx/gi;

#my $workbook = Spreadsheet::WriteExcel->new($xls_file);
my $workbook = Excel::Writer::XLSX->new($xls_file);
my $worksheet = $workbook->add_worksheet("human gene families");
open(FH,"< $hgnc_family_file") or die "Cannot open file: $!\n";
my ($x,$y) = (0,0);
while (<FH>){
 chomp;
 my @list = split /\t/,$_;
 foreach my $c (@list){
    $worksheet->write($x, $y++, $c);
 }
 $x++;$y=0;
}
close(FH);
$workbook->close();

exit 0;
