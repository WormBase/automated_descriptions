#!/usr/bin/env perl
use ConciseDescriptions;
use File::Slurp; 
use List::MoreUtils qw(uniq);
use LWP::Simple;
use LWP::UserAgent;
use strict;
use warnings;
# 
my $web_page = "http://tazendra.caltech.edu/~azurebrd/cgi-bin/forms/generic.cgi?action=WpaXrefBackwards";
my $wbpaperhash_ref=fillwbpaperhash($web_page);
my %wbpaperhash = %$wbpaperhash_ref;
my $double_space = "  ";
my $space = " ";
my $dash_line="";
# html path is defined by reading the an input file (e.g., docroot in RH, /var/www in Ubuntu, etc.)
my $html = ConciseDescriptions::get_html_dir();
my $RELEASE = ConciseDescriptions::get_release();
my $PRODUCTION_RELEASE =  ConciseDescriptions::get_production_release();
#
my $species = $ARGV[0];
my $project = $ARGV[1];
my $species_name = $ARGV[2];
my $species_prefix = $ARGV[3];

chomp($species);
chomp($project);
chomp($species_name);

if ($species_prefix){
 chomp($species_prefix);
 $species_prefix =~ s/^\s+//;
 $species_prefix =~ s/\s+$//;
}

$project =~ s/^\s+//;
$project =~ s/\s+$//;
$species =~ s/^\s+//;
$species =~ s/\s+$//;
$species_name =~ s/^\s+//;
$species_name =~ s/\s+$//;
#
my $home = $html . "concise_descriptions/";
my $path = $home . "release/$PRODUCTION_RELEASE/$species/descriptions/";
my $description_directory= $path . "individual_gene_descriptions/";
my $orthology = $home . "release/$PRODUCTION_RELEASE/$species/orthology/";
my $elegans_orthology = $home . "release/$PRODUCTION_RELEASE/c_elegans/orthology/";
my $summary = $path . "OA_concise_descriptions.txt";
# Create advertising string
my $advertising_string = "This description was generated automatically by a Textpresso script based on homology/orthology data and Gene Ontology (GO) annotations from the $RELEASE version of WormBase.";
# Remove previous incarnation of summary file
if (-e $summary){
  my @args = ("rm", "-f", $summary);
  system(@args) == 0 or die("could not delete file $summary\n");
}
#
my $current_date = `date +"%Y-%m-%d"`;
$current_date = substr($current_date,0,-1);
#
# BioMart file
#
my $biomart = $elegans_orthology . "input_files/HumanIDs_mart_export.txt";
my $reading_biomart_file_for_error = read_file($biomart);
if ($reading_biomart_file_for_error =~ /ERROR/){
    $biomart =  $elegans_orthology . "input_files/HumanIDs_mart_export.org.txt";
}
#
my $xrefs = $orthology . "input_files/$species.$RELEASE.xrefs.txt";
my $blastp = $orthology . "input_files/$species.$RELEASE.best_blastp_hits.txt";
my $gene_protein_hash_ref = ConciseDescriptions::get_gene_protein_hash($species, $xrefs); 
my %gene_protein_hash = %$gene_protein_hash_ref;
my $species_ensp_hash_ref = ConciseDescriptions::get_species_ensp_hash($species, $blastp);
my %species_ensp_hash = %$species_ensp_hash_ref;
my $gene_association_file ="ftp://ftp.sanger.ac.uk/pub/wormbase/releases/$RELEASE/ONTOLOGY/gene_association.$RELEASE.wb.$species";
my $gene_association_contents = ConciseDescriptions::getwebpage($gene_association_file);
my @gene_association_lines = split(/\n/, $gene_association_contents);
    my %reference=();
    my %accession=();
    foreach my $line2 (@gene_association_lines) {
      my $gene2 = $line2;
      chomp($gene2);
      next if ($gene2 =~ m/\!/);
      my @gene_array2 = split(/\t/,$gene2);     
      my $gene_id2 = $gene_array2[1];
      $reference{$gene_id2}="";
      $accession{$gene_id2}="";
    }
    foreach my $line (@gene_association_lines) {
      my $gene = $line;
      chomp($gene);
      next if ($gene =~ m/\!/);
      
      my @gene_array = split(/\t/,$gene);
      
      my $gene_id = $gene_array[1];
      my $gene_name = $gene_array[2];
      my $GOID   = $gene_array[4];
      my $ref = $gene_array[5];
      my $evidence_code = $gene_array[6];
      my $with = $gene_array[7];

      chomp($gene_id);
      chomp($gene_name);
      chomp($GOID);
      chomp($evidence_code);
      chomp($ref);
      chomp($with);
 
     $with =~s/^\s+//;
     $with =~s/\s+$//;
     $ref =~s/^\s+//;
     $ref =~s/\s+$//;
     $ref =~s/\t+$//;
     $ref =~s/^\t+//;
     $with =~s/\t+$//;
     $with =~s/^\t+//;

      if (($with !~m/WB\:WBRNAi/) or ($with !~ m/WBPhenotype/)) {
      if ($accession{$gene_id} ne "") {
       if ($accession{$gene_id} !~ m/$with/){
          $accession{$gene_id} .= "\," . $with;
       }
      } else {
          $accession{$gene_id} = $with;
      }
      }
      if ($reference{$gene_id} ne ""){
        if ($reference{$gene_id} !~ m/$ref/){
          $reference{$gene_id} .= "\|" . $ref;
        }
      } else {
          $reference{$gene_id} = $ref;
      }
} # end foreach @gene_association_lines

#
my @description_files  = glob("$description_directory/WBGene*");
my $description="";
#
# Create a list of genes
#
my @list = (@description_files);
#
 my @unsorted_files = ();
 foreach my $item (@list){
     my @items = split("/",$item);
     my $size_of_array = @items;
     my $file = $items[$size_of_array-1];
     push(@unsorted_files, $file);
#    $item =~s/^[^WBG]*(?=WBG)//;    
#    print "$item\n";
 }
 my @sort_unique_list = uniq(sort(@unsorted_files));

foreach my $file (@sort_unique_list){
 next if ($file eq "WBGene00000000");
 my $gene_id = $file;
 my $description_file = $description_directory . $file;

 if (-e $description_file) {
     $description = read_file($description_file);
                        }
                        
     $description =~s/^\s+//;
     $description =~s/\s+$//;

 my $ref="";
 my $acc="";
 if ($reference{$gene_id}){
     $ref = $reference{$gene_id};
 }

 if ($accession{$gene_id}){
     $acc = $accession{$gene_id};
 }
     $acc =~s/^\s+//;
     $acc =~s/\s+$//;
     $ref =~s/^\s+//;
     $ref =~s/\s+$//;

     $ref =~s/WB\_REF\://g;

     $ref=~s/^\,+//g;
     $acc=~s/^\,+//g;
     $ref=~s/\,+$//g;
     $acc=~s/\,+$//g;
     $ref=~s/\t+/\,/g;
     $ref=~s/\s+/\,/g;
     $ref=~s/\|/\,/g;     
     
     my $ensemble = "";
     my $ensp = "";
     my $gene_protein = $gene_protein_hash{$gene_id};
     if ($gene_protein){
        $ensp = $species_ensp_hash{$gene_protein};
     }
     if ($ensp){
        $ensemble = $ensp;
     }

# $acc =~s/\,/ \,/g;

 my @refs = split(/\,/, $ref);
 my @new_refs=();
 foreach my $r (@refs){
  if ($r =~/PMID\:/){
      $r=~/PMID\:(\d+)/;
      my $pmid = $1;
      my $wb = $wbpaperhash{$pmid};
      if ($wb){
          $r = $wb;
      } else {
      print "No WBPaperID for\tPMID\:$pmid\n";
    }
  }
  push(@new_refs, $r);
 }
 my @uniq_refs = uniq(sort(@new_refs)); 
 $ref = join("\,",@uniq_refs);
 $ref =~ s/\,(?=\S)/\, /g;
 $acc =~ s/\,(?=\S)/\, /g;
 $ref =~ s/GO_REF:0000002/WBPaper00045688, WBPaper00045689/g;
 $acc =~s/\|/\, /g;
 my $output;

 if (($acc) and ($ensemble)) {
 $output = $gene_id . "\t" . $current_date . "\t" . $ref . "\t" . $ensemble . "\, " . $acc . "\t" . $description;} elsif ($ensemble) {
 $output = $gene_id . "\t" . $current_date . "\t" . $ref . "\t" . $ensemble . "\t" . $description;
} else {
$output = $gene_id . "\t" . $current_date . "\t" . $ref . "\t" . "\t" . $description;
}

    $output =~ s/  / /g;
    $output =~ s/$double_space/$space/g;

    $output =~ s/^\s+//;
    $output =~ s/\s+$//;

my $last = chop($output);
if ($last=~/\;/){
    $output .= "\.\t$species_name\t$advertising_string\n";
} else{
    $output =~ s/\n$//;
    $output =~ s/\;$//;
    $output .= "\.\t$species_name\t$advertising_string\n";
}
 
 write_file($summary, {append => 1 }, $output);

}
#
exit 0;

sub fillwbpaperhash{
    my $web_page = shift;
    my $contents = ConciseDescriptions::getwebpage($web_page);
    my @lines = split(/\n/, $contents);
my %wbpaperhash=();
    for (my $i=0; $i<@lines; $i++) {
        if ($lines[$i] =~ /pmid/) {
            $lines[$i] =~ /pmid(\d+)/;
            my $pmid = $1;
            $lines[$i] =~ /WBPaper(\d+)/;
            my $wbpaper_id = "WBPaper" . $1;
            if (($wbpaper_id) and ($pmid)) {
                $wbpaperhash{$pmid}=$wbpaper_id;
            }
        }
    }
return \%wbpaperhash;
}
