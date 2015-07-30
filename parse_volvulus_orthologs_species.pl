#!/usr/bin/env perl
use File::Slurp; 
use strict;
use warnings;
use Data::Dumper;
use ConciseDescriptions;
#
# Creates a homolog file for WBGenes, named $output_file.
# Authors: J. Done and R. Kishore, California Institute of Technology, 2014. 
#
# $file should be pointing to the latest ortholog file
#
#
my $species = $ARGV[0];
chomp($species);

$species =~ s/^\s+//;
$species =~ s/\s+$//;

if (($species =~/malayi/) or ($species =~/ratti/)) { 

my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";

my $orthology = $home . "release/$PRODUCTION_RELEASE/$species/orthology/input_files/";

my $output_file = $orthology . "orthologs.$species.volvulus.txt";
my $file = $orthology . "$species.orthologs.txt";
# if the output file exists delete it.
if (-e $output_file){
   my @args = ("/bin/rm", "-f", $output_file);
   system(@args) == 0 or die("could not delete file $output_file\n");
}
my @ortholog_lines = read_file($file); 
my $wbline;
my $ortholog_line;
my $oldgene="";
foreach my $line (@ortholog_lines){
  if ($line =~/^WBGene/){
   chomp($line);
   $wbline = $line;
   $ortholog_line="";
  }
  if ($line =~/volvulus/){
   chomp($line);
   $ortholog_line = $wbline . " " . $line;
   my $newgene = $wbline;
   if ($newgene ne $oldgene){
   write_file($output_file, {append => 1 }, "\n");
   write_file($output_file, {append => 1 }, "\n");
   $oldgene=$newgene;
  }
   write_file($output_file, {append => 1 }, $ortholog_line);
   write_file($output_file, {append => 1 }, "\n");
  } 

}

}
exit 0;
