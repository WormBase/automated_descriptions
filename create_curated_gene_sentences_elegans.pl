#!/usr/bin/env perl
use warnings;
use strict;
use ConciseDescriptions;
use File::Slurp; 
use List::MoreUtils qw(uniq);

my $species="c\_elegans";
my $species_name = "Caenorhabditis elegans";

my $html = ConciseDescriptions::get_html_dir();
my $RELEASE = ConciseDescriptions::get_release();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $concise_descriptions = $html . "concise_descriptions/";
my $individual_path = $concise_descriptions . "release/$PRODUCTION_RELEASE/$species/descriptions/individual_gene_descriptions/";
my $output_path = $concise_descriptions . "release/$PRODUCTION_RELEASE/$species/descriptions/";
my $summary = $output_path . "WBGenes_concise_descriptions_for_manual.txt";
my $elegans_gene_list_dir = $concise_descriptions . "release/$PRODUCTION_RELEASE/c_elegans/gene_lists/";
# if the output file exists delete it.
if (-e $summary){
   my @args = ("/bin/rm", "-f", $summary);
   system(@args) == 0 or die("could not delete file $summary\n");
}
   write_file($summary,"\n");

#
my $curated_gene_list = $elegans_gene_list_dir . "sort.curated_genes.txt";
my @curated_genes = read_file($curated_gene_list);

my @individual_files = glob("$individual_path/WBGene*");
   foreach my $file (@individual_files){
    if (-e $file){
    my $gene_id = $file;
       chomp($gene_id);
       $gene_id =~ s{\.[^.]+$}{};
       $gene_id =~ s{.*/}{};
       $gene_id =~ s/^\s+//;
       $gene_id =~ s/\s+$//; 
    if (grep {$_ =~/$gene_id/} @curated_genes) {
      print "...processing $gene_id\n";
      my $sentence = read_file($file);
      write_file($summary, {append => 1 }, $gene_id);
      write_file($summary, {append => 1 }, "\n");
      write_file($summary, {append => 1 }, $sentence);
      write_file($summary, {append => 1 }, "\n\n");
    }
   }
}
exit 0;
