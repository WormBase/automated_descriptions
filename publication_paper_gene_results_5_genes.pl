#!/usr/bin/env perl
use strict;
use warnings;
use ConciseDescriptions;
use List::MoreUtils qw(uniq);
use File::Slurp;
use LWP::Simple;

my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/c_elegans/";
my $gene_list_dir = $home . "gene_lists/";
my $input_dir = $html . "concise_descriptions/textpresso/";
my $output_dir = $html . "concise_descriptions/textpresso/";
my $input = $input_dir . "textpresso_papers_results_genes.all.txt";
my $output = $output_dir . "textpresso_papers_results_genes.txt";
if (-e $output){
  my @args = ("rm", "-f", $output);
  system(@args) == 0 or die("could not delete file $output\n");
}
my $header="  paper id\t\tgenes\n";
write_file($output, {append => 1 }, $header);
 
 if (-e $input){

 my @lines = read_file($input);
  foreach my $line (@lines){
  next if ($line=~/paper id/);
  chomp($line);
  my $paper = "";
  my $list = "";
  $line =~/(WBPaper\S+)/;
  $paper = $1;
  print "$paper\n";
  $line =~s/$paper//g;
  $list  = $line;
     $list =~s/^\s+//;
     $list =~s/\s+$//;
  my @genes = split(/\s/, $list);
  if (@genes < 6){ 
   my $output_string = $paper . " " . $list . "\n";   
   write_file($output, {append => 1 }, $output_string);
  }
  }
}

exit 0;
