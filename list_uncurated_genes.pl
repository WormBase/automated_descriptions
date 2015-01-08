#!/usr/bin/env perl
use warnings;
use strict;
use File::Slurp;
use List::MoreUtils qw(uniq);
use ConciseDescriptions;

my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
my $go_elegans_dir = $home . "release/$PRODUCTION_RELEASE/c_elegans/gene_ontology/input_files/";
my $gene_association_file = $go_elegans_dir . "gene_association.$RELEASE.wb.c_elegans";
my $wbgene_gaf_elegans_array_ref = get_wbgene_gaf_array($gene_association_file);
my @wbgene_gaf_elegans_array = @$wbgene_gaf_elegans_array_ref;
my $gene_list_dir = $home . "release/$PRODUCTION_RELEASE/c_elegans/gene_lists/";
my $curated_gene_list = $gene_list_dir . "sort.curated_genes.txt";
my @curated_genes = read_file($curated_gene_list);
my $elegans_orthology = $home . "release/$PRODUCTION_RELEASE/c_elegans/orthology/input_files/";
my $gene_orthologs_file = $elegans_orthology . "c_elegans.PRJNA13758.$RELEASE.orthologs.txt";
my $wbgene_ortholog_elegans_array_ref = get_wbgene_ortholog_array($gene_orthologs_file);
my @wbgene_ortholog_elegans_array = @$wbgene_ortholog_elegans_array_ref;

my $output_file =  $gene_list_dir . "sort.uncurated_genes.txt";
# if the output file exists delete it.
if (-e $output_file){
   my @args = ("/bin/rm", "-f", $output_file);
   system(@args) == 0 or die("could not delete file $output_file\n");
}
# combine gene lists
my @gene_array = @wbgene_ortholog_elegans_array;
foreach my $element (@wbgene_gaf_elegans_array){
 push(@gene_array, $element);
}
my @sorted_unique_gene_array = sort(uniq(@gene_array));

# subtract curated genes
my @uncurated_genes_array = ();
foreach my $test (@sorted_unique_gene_array){
 my $keep =0;
 foreach my $curated (@curated_genes){
   if ($curated =~/$test/){
    $keep = 1;
   }
 }
 if ($keep ==0){
   push(@uncurated_genes_array, $test);
 }
}
# print remaining genes
my @sorted_uncurated_genes_array = sort(@uncurated_genes_array);

foreach my $element (@sorted_uncurated_genes_array){
        my $out = $element . "\n";
        write_file($output_file, {append => 1 }, $out);
} 

exit 0;

sub get_wbgene_gaf_array{
 my $file = shift;
 my @array=();
 my @lines = read_file($file);
  foreach my $line (@lines){
  my @fields = split(/\t/, $line);
  my $gene_name = $fields[2];
  my $gene_id = $fields[1];
  if ($gene_name){
       push(@array, $gene_id);
   }
 }
 return \@array;
}
sub get_wbgene_ortholog_array{
 my $file = shift;
 my @array=();
 my @lines = read_file($file);
  foreach my $line (@lines){
   next if ($line =~/\#/);
   next if ($line =~/\=/);
   
   my @fields = split(/\t/, $line);
   my $size = @fields;
   next if ($size gt 2);
#   print "size is $size\n";
#   if ($size==2){
#       print "$fields[0]\t$fields[1]\n";
#   }
   my $gene_name = $fields[1];
   my $gene_id = $fields[0];
   if ($gene_id=~/WBGene/){
       push(@array, $gene_id);
   }
  }
 return \@array;
}
