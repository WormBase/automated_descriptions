#!/usr/bin/env perl 
use Carp;
use LWP::Simple;
use File::Slurp;
use ConciseDescriptions;
use strict;
use warnings;
#
# Authors: J. Done and R. Kishore, California Institute of Technology, 2014. 
# http://wiki.wormbase.org/index.php/Generation_of_automated_descriptions
#
# html path is defined by reading the an input file (e.g., docroot in RH, /var/www in Ubuntu, etc.)
my @args;
my $status;
my $AND = "AND";
my $species_project = $ARGV[0];
 chomp($species_project);
my ($species, $prefix) = split(/$AND/, $species_project);
 $species =~s/^\s+//;
 $species =~s/\s+$//;
 if ($prefix){
  $prefix =~s/^\s+//;
  $prefix =~s/\s+$//;
  $prefix = lc $prefix;
 }
 if ($species =~/elegans/){
  $prefix = "ce";
 }
#
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
my $RELEASE = ConciseDescriptions::get_release();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
#
my $anatomy_ec_dir = $home . "release/$PRODUCTION_RELEASE/$species/anatomy\_expression\_cluster/input\_files/";
my $mole_reg_ec_dir = $home . "release/$PRODUCTION_RELEASE/$species/molecule\_regulation\_expression\_cluster/input\_files/";
my $gene_reg_ec_dir = $home . "release/$PRODUCTION_RELEASE/$species/gene\_regulation\_expression\_cluster/input\_files/";

my $anatomy_file  = $anatomy_ec_dir . $prefix . "ECsummary_anatomy.$RELEASE.txt";
my $molecule_file = $mole_reg_ec_dir . $prefix . "ECsummary_molReg.$RELEASE.txt";
my $gene_file = $gene_reg_ec_dir . $prefix . "ECsummary_geneReg.$RELEASE.txt";

my $url = "ftp://caltech.wormbase.org/pub/wormbase/ExprClusterSummary/$PRODUCTION_RELEASE/";
my $url_anatomy = $url . $prefix . "ECsummary_anatomy.$RELEASE.txt";
my $url_molecule = $url . $prefix . "ECsummary_molReg.$RELEASE.txt";
my $url_gene = $url . $prefix . "ECsummary_geneReg.$RELEASE.txt";
#
if (head($url_anatomy)) {
    getstore($url_anatomy, $anatomy_file);
}
#
if (head($url_molecule)) {
    getstore($url_molecule, $molecule_file);
}
#
if (head($url_gene)) {
    getstore($url_gene, $gene_file);
}
#
exit 0;
