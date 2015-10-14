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
my $AND = "AND";
my $species_term = $ARGV[0];
my ($species, $project, $species_name, $species_prefix) = split(/$AND/, $species_term);

chomp($species);
chomp($project);
chomp($species_name);

 $species_name =~s/\_/ /g;

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
my $go_home = $home . "release/$PRODUCTION_RELEASE/$species/gene_ontology/";
my $tissue_home = $home . "release/$PRODUCTION_RELEASE/$species/tissue_expression/";
my $summary = $path . "OA_concise_descriptions.$PRODUCTION_RELEASE.txt";
my $number_file = $path . "OA_concise_descriptions.number_of_genes.$PRODUCTION_RELEASE.txt";
# Create advertising string
my $advertising_string = "This description was generated automatically by a Textpresso script based on homology/orthology data, Gene Ontology (GO) annotations and tissue expression data from the $RELEASE version of WormBase.";
# Remove previous incarnation of summary and number files
if (-e $summary){
  my @args = ("rm", "-f", $summary);
  system(@args) == 0 or die("could not delete file $summary\n");
}
my $count = 0;
if (-e $number_file){
  my @args = ("rm", "-f", $number_file);
  system(@args) == 0 or die("could not delete file $number_file\n");
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

my $elegans_sapiens_orthology = $elegans_orthology . "input_files/orthologs.Caenorhabditis_elegans.H_Sapiens.txt";
my ($gene_ensg_hash_ref, $gene_hgnc_hash_ref) = get_sapiens_hash($elegans_sapiens_orthology);
my %gene_ensg_hash= %{$gene_ensg_hash_ref};
my %gene_hgnc_hash= %{$gene_hgnc_hash_ref};

my $gene_association_file = "";
my @tissue_gene_association_lines = ();
my $tissue_gene_association_file = "";
if ($species=~/elegans/){
    $gene_association_file = $go_home . "input_files/gene_association.wb";
    $tissue_gene_association_file = $tissue_home . "input_files/anatomy_association.$RELEASE.wb";
    @tissue_gene_association_lines = read_file($tissue_gene_association_file);
} else{
    $gene_association_file = $go_home . "input_files/gene_association.$RELEASE.wb.$species";
}
 my @gene_association_lines = read_file($gene_association_file);
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

#     print "$gene_id\tref with\t$ref\t$with\n";
 
     $with =~s/^\s+//;
     $with =~s/\s+$//;
     $ref =~s/^\s+//;
     $ref =~s/\s+$//;
     $ref =~s/\t+$//;
     $ref =~s/^\t+//;
     $with =~s/\t+$//;
     $with =~s/^\t+//;

     my @with_array = split(/\|/, $with);
     my $new_with = "";
     foreach my $w (@with_array){
      if ($w =~/InterPro/){
       $new_with .= $w ."\|";
      }
     }
       $new_with =~s/\|+$//g;
       $with = $new_with;
#     print "with is $with\n";
  
     $ref =~ s/GO_REF:0000002/WBPaper00045688|WBPaper00045689/g; 
     $ref =~s/WB\_REF\://g;

#     print "ref is $ref\n";
     if (($with =~ /InterPro/) or ($with =~/PMID/i) or ($with =~/WBPaper/))  {
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

    foreach my $line (@tissue_gene_association_lines) {
      my $gene = $line;
      chomp($gene);
      next if ($gene =~ m/\!/);
      
      my @gene_array = split(/\t/,$gene);
      
      my $gene_id = $gene_array[1];
      my $gene_name = $gene_array[2];
      my $ref = $gene_array[5];

      chomp($gene_id);
      chomp($gene_name);
      chomp($ref);

     $ref =~s/^\s+//;
     $ref =~s/\s+$//;
     $ref =~s/\t+$//;
     $ref =~s/^\t+//;

     $ref =~ s/GO_REF:0000002/WBPaper00045688|WBPaper00045689/g; 
     $ref =~s/WB\_REF\://g;

<<<<<<< HEAD
      if (exists $reference{$gene_id}){
=======
      if ($reference{$gene_id} ne ""){
>>>>>>> b8516fd33e3bbe9e4fd676bcfd68cf86d33ceb53
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

  if ($species =~/elegans/){
  my $ensg_list = "";
  if ($gene_ensg_hash{$gene_id}){
      $ensg_list = "ENSEMBL\:" . $gene_ensg_hash{$gene_id};
   if ($ensg_list =~/AND/){
      $ensg_list =~s/AND/\, ENSEMBL\:/g;
    }
      $ensemble =  $ensg_list;
    }
   }

# $acc =~s/\,/ \,/g;

 my @refs = split(/\,/, $ref);
 my @new_refs=();
 foreach my $r (@refs){
  if (($r =~/GO/) or ($r=~/WBPaper/) or ($r =~/PMID/)){
  $r =~ s/GO_REF:0000002/WBPaper00045688, WBPaper00045689/g;  
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
  if ($r =~/WBPaper/i){
   push(@new_refs, $r);
   }
  }
 }
 my @uniq_refs = uniq(sort(@new_refs)); 
 $ref = join("\,",@uniq_refs);
 $ref =~ s/\,(?=\S)/\, /g;
 $acc =~ s/\,(?=\S)/\, /g;
 $acc =~s/\|/\, /g;
 my $output;
 
if (($acc) and ($ensemble)) {
 $output = $gene_id . "\t" . $current_date . "\t" . $ref . "\t" . $ensemble . "\, " . $acc . "\t" . $description;} 
 elsif ($ensemble) {
 $output = $gene_id . "\t" . $current_date . "\t" . $ref . "\t" . $ensemble . "\t" . $description;
 } elsif  ($acc) {
  $output = $gene_id . "\t" . $current_date . "\t" . $ref . "\t" . $acc . "\t" . $description;} 
  else {
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

 $count++; 
 write_file($summary, {append => 1 }, $output);

}
#
 my $number_of_genes = "The number of genes written in the OA tab delimited file for $species_name is $count\.\n";
 write_file($number_file, $number_of_genes);
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
sub get_sapiens_hash{
 my $orthology=shift;
 my @lines = read_file($orthology);
 my $old_gene="";
 my %gene_ensg_hash=();
 my %gene_hgnc_hash=();
 my @multiple_source_array=();
 my @single_source_array=();
 my $oldgene="";
 my $newgene; 
 foreach my $line (@lines){    
    if ($line =~/WBGene/){
     my ($wb_gene_id, $gene_name, $species_name, $ensg, $hgnc, $source) = split(/\t/, $line); 
        $wb_gene_id =~ s/^\s+//;
        $wb_gene_id =~ s/\s+$//;
        $newgene = $wb_gene_id;
     if ($newgene ne $oldgene){
         my $single_size = @single_source_array;
         my $multiple_size = @multiple_source_array;
#         print "multiple size is $multiple_size\n";
#         print "single size is $single_size\n";

         if ($multiple_size gt 0){
             foreach my $a_line (@multiple_source_array){
#	         if ($gene_name =~/F58E6\.13/){
#                  print " a line is $a_line\n";
#         	}       
               my ($a_wb_gene_id, $a_gene_species_name, $a_ensg, $a_hgnc, $a_source) = split(/\t/, $a_line);
               my ($a_gene_name, $a_species_name) = split(/ /, $a_gene_species_name);
                   $gene_hgnc_hash{$a_ensg}  = $a_hgnc;

               if ($gene_ensg_hash{$a_wb_gene_id}){
                   $gene_ensg_hash{$a_wb_gene_id} .= "AND" .  $a_ensg;
               } else{
                   $gene_ensg_hash{$a_wb_gene_id} = $a_ensg;
               }
            }
          } else {
             foreach my $a_line (@single_source_array){
#	         if ($gene_name =~/F58E6\.13/){
#                  print " a line is $a_line\n";
#         	}       
               my ($a_wb_gene_id, $a_gene_species_name, $a_ensg, $a_hgnc, $a_source) = split(/\t/, $a_line);
               my ($a_gene_name, $a_species_name) = split(/ /, $a_gene_species_name);
                   $gene_hgnc_hash{$a_ensg}  = $a_hgnc;

               if ($gene_ensg_hash{$a_wb_gene_id}){
                   $gene_ensg_hash{$a_wb_gene_id} .= "AND" .  $a_ensg;
               } else{
                   $gene_ensg_hash{$a_wb_gene_id} = $a_ensg;
               }
          }
         }
         @multiple_source_array=();
         @single_source_array=();
         $oldgene=$newgene;
     }
         
         if ($line =~/\;/){
          my $multiple = @multiple_source_array;
          if ($multiple lt 5){
           push(@multiple_source_array, $line);
          }
       } else {
          my $single = @single_source_array;
          if ($single lt 5){
           push(@single_source_array, $line);
          }
       }
  }
 }
 return \%gene_ensg_hash, \%gene_hgnc_hash;
}
