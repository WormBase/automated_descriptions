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
#my $output_dir = "./";
#
my $dir = $html . "concise_descriptions/release/$RELEASE/c_elegans/";

my $deadList = $dir . "gene_lists/sort.dead_genes.txt";
my $geneIDs =  $dir . "gene_lists/geneIDs.txt";

my @descriptions = glob( $dir . 'descriptions/individual_gene_descriptions/*' );
my @molecular_functions = glob( $dir . 'gene_ontology/output_files/individual_function_sentences/*' );
my @components = glob( $dir . 'gene_ontology/output_files/individual_component_sentences/*' );
my @processes = glob( $dir . 'gene_ontology/output_files/individual_process_sentences/*' );
my @dead_genes = read_file($deadList); 

#print "descriptions \= @descriptions\n";

my %dead_hash = ();
my %process_hash = ();
my %function_hash = ();
my %component_hash = ();

foreach my $d (@dead_genes) {
 if ($d){
     # remove leading whitespace
     $d =~ s/^\s+//;
     # remove trailing whitespace
     $d =~ s/\s+$//;
  $dead_hash{$d} = $d;
 }
}
foreach my $p (@processes) {
 if ($p){
     $p =~ s/^.*[\/\\]//;
     # remove leading whitespace
     $p =~ s/^\s+//;
     # remove trailing whitespace
     $p =~ s/\s+$//;
#  print "$p\n";
  $process_hash{$p} = $p;
 }
}
foreach my $f (@molecular_functions) {
 if ($f){
     $f =~ s/^.*[\/\\]//;
     # remove leading whitespace
     $f =~ s/^\s+//;
     # remove trailing whitespace
     $f =~ s/\s+$//;
  $function_hash{$f} = $f;
 }
}
foreach my $c (@components) {
 if ($c){
     $c =~ s/^.*[\/\\]//;
     # remove leading whitespace
     $c =~ s/^\s+//;
     # remove trailing whitespace
     $c =~ s/\s+$//;
  $component_hash{$c} = $c;
 }
}
my %gene_ids = ();
my @genes = ();
my %gene_string = ();
my %gene_names = ();
my @lines = read_file($geneIDs);
foreach my $line (@lines){
 chomp($line);
 next if ($line =~/dead/i);
 my ($no, $gene_id, $gene_name, $alt_name, $live) = split(/\,/, $line);
 next if ($component_hash{$gene_id});
 next if ($process_hash{$gene_id});
 next if ($function_hash{$gene_id});
 next if ($dead_hash{$gene_id});

 if ($gene_id){
  push(@genes, $gene_id);
  if ($gene_name){
    $gene_ids{$gene_name} = $gene_id;
    $gene_names{$gene_id} = $gene_name;
    $gene_string{$gene_name} = $gene_name;
  } elsif ($alt_name){
    $gene_ids{$alt_name} = $gene_id;
    $gene_names{$gene_id} = $alt_name;
    $gene_string{$alt_name} = $alt_name;
  }
 }
}
#
my %hash_pubs = ();
my %hash_string = ();
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
#  print "$paper\n";
  $line =~s/$paper//g;
  $list  = $line;
     $list =~s/^\s+//;
     $list =~s/\s+$//;
  my $newlist = "";
  my @genes = split(/\s/, $list);
  my $count = 0;
  if (@genes < 6){
   foreach my $gene (@genes){
    $gene =~/(.*)\(/;
    my $name = $1;
    if ($gene_string{$name}) {
#    print "$name\n";
    $gene =~/\S+\((.*)\)/;
    my $number = $1;
#    print "number is $number\n";
    $count+= $number;
    $newlist .= " " . $name . "\(" . $number . "\)"; 
    }
   }
     $newlist =~s/^\s+//;
     $newlist =~s/\s+$//;
   my $output_string = "";
#   print "score is $count\n";
   if ($newlist){
       $output_string = $paper . " " . $newlist . "\n";
       $hash_pubs{$paper} = $count;
       print "$paper\t$count\n";
       $hash_string{$paper} = $output_string; 
   }  
#   write_file($output, {append => 1 }, $output_string);
  }
  }

foreach (reverse sort { ($hash_pubs{$a} <=> $hash_pubs{$b}) || ($a cmp $b) } keys %hash_pubs) 
{
    my $key = $_;
    my $out = $hash_string{$key};
    write_file($output, {append => 1 }, $out);
#    print "$_: $userids{$_}\n";
}
}
exit 0;
