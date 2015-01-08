#!/usr/bin/env perl
use ConciseDescriptions;
use List::MoreUtils qw(uniq);
use List::Util qw/first/;
use File::Slurp;
use Text::CSV;
use warnings;
use strict;

my $species = $ARGV[0];
my $project = $ARGV[1];

chomp($species);
chomp($project);
           
$project =~ s/^\s+//;
$project =~ s/\s+$//;
$species =~ s/^\s+//;
$species =~ s/\s+$//;

my $RELEASE = ConciseDescriptions::get_release();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";

my $orthology = $home . "release/$PRODUCTION_RELEASE/$species/orthology/";
my $elegans_orthology = $home . "release/$PRODUCTION_RELEASE/c_elegans/orthology/";
my $output_file = $orthology . "output_files/sentences_for_orthology.txt";
#
if (-e $orthology){
   print "$orthology exists\n"; 
} else {
   mkdir $orthology or die "could not create $orthology";
}
my $individual_path = $orthology . "output_files/individual_gene_sentences/";
if (-e $individual_path){
   print "$individual_path exists\n"; 
} else {
   mkdir $individual_path or die "could not create $individual_path";
}
# if the output file exists delete it.
if (-e $output_file){
   my @args = ("/bin/rm", "-f", $output_file);
   system(@args) == 0 or die("could not delete file $output_file\n");
}
# if the output files for individual gene sentences exist delete them.
   my @individual_files = glob("$individual_path/WBGene*");
      foreach my $individual (@individual_files){
      my @args = ("/bin/rm", "$individual");
      system(@args) == 0 or die("could not delete $individual\n");
   }

my $xrefs = $orthology . "input_files/$species.$RELEASE.xrefs.txt";
my $blastp = $orthology . "input_files/$species.$RELEASE.best_blastp_hits.txt";

my $biomart = $elegans_orthology . "input_files/HumanIDs_mart_export.txt";
my $reading_biomart_file_for_error = read_file($biomart);
if ($reading_biomart_file_for_error =~ /ERROR/){
    $biomart =  $elegans_orthology . "input_files/HumanIDs_mart_export.org.txt";
}
my $xrefs_elegans = $elegans_orthology . "input_files/c_elegans.$RELEASE.xrefs.txt";
my $orthology_elegans = $elegans_orthology . "input_files/c_elegans.$project.$RELEASE.orthologs.txt";
my $wbgene_id_human = $orthology . "input_files/wb_gene_id_human.txt";

my $gene_name_hash_ref = get_gene_name_hash($species, $xrefs);
my %gene_name_hash = %$gene_name_hash_ref;
my $gene_protein_hash_ref = ConciseDescriptions::get_gene_protein_hash($species, $xrefs); 
my %gene_protein_hash = %$gene_protein_hash_ref;
my $species_ensp_hash_ref = ConciseDescriptions::get_species_ensp_hash($species, $blastp);
my %species_ensp_hash = %$species_ensp_hash_ref;
my $species_ce_protein_hash_ref = get_species_ce_protein_blastp_hash($species, $blastp);
my %species_ce_protein_hash =%$species_ce_protein_hash_ref;
my $elegans_gene_name_hash_ref = get_gene_name_hash("c_elegans", $xrefs_elegans);
my %elegans_gene_name_hash = %$elegans_gene_name_hash_ref;
my $elegans_gene_protein_hash_ref = get_gene_protein_hash("c_elegans", $xrefs_elegans);
my %elegans_gene_protein_hash = %$elegans_gene_protein_hash_ref;
my $biomart_description_hash_ref = get_biomart_description_hash($biomart);
my %biomart_description_hash = %$biomart_description_hash_ref;
my $biomart_hgnc_hash_ref = get_biomart_hgnc_hash($biomart);
my %biomart_hgnc_hash = %$biomart_hgnc_hash_ref;
my $gene_array_ref=get_gene_array($xrefs);
my @gene_array = @$gene_array_ref;
my $elegans_name_hash_ref = get_inverse_gene_protein_hash("c_elegans", $xrefs_elegans);
my %elegans_name_hash = %$elegans_name_hash_ref;

my $elegans_string = " encodes an ortholog of C. elegans ";
my $elegans_string1 = " encodes an ortholog of C. elegans ";
my $and = "and";
my $comma = "\,";
my $and_of = " and an ortholog of ";
my $and_human_string = " and human protein ";
my $only_human_string = " encodes an ortholog of human ";
my $only_human_string1 = " encodes an ortholog of ";
my $family = "\, family ";
my $family2 = " family member";
my $family3 = " family";
my $family4 = " family\, member";
my $uncharacterized = "Uncharacterized";
my $group = "\, group ";
my $member = "\, member ";
my $subfamily = "\, subfamily ";
my $subfamily2 = "\, sub\-family ";
my $subfamily3 = "subfamily";
my $polypeptide ="\, polypeptide ";
my $class = "\, class ";
my $doublespace = "  ";
my $space = " ";
my $spacesemicolon = " \;";
my $semicolon = "\;";
my $homolog  = "homolog ";
my $homolog2 = "homolog";
my $homolog3 = "homolog\.*";

# Define the various model organisms used
my @mods=("S\. cerevisiae","yeast","S\. pombe","mouse","Drosophila","Chlamydomonas","zebrafish","Rat","rat","chicken","Strongylocentrotus purpuratus","Arabidopsis","Dictyostelium","bacterial","E\. coli","H\. influenzae","Dictyostelium");
 my $elegans = "C\. elegans";
 my $human = " and human ";
 my $sentence  = "";
 my $ensp = "";
 my $hgnc = "";
 my $elegans_protein_id ="";
 my $description = "";
 my $elegans_ortholog_gene = "";
foreach my $gene_id (uniq @gene_array){
    $description ="";
    $elegans_protein_id = "";
    $elegans_ortholog_gene = "";
# print "$gene_id\n";
 my $out = $individual_path . $gene_id;
 my $gene_name = $gene_name_hash{$gene_id};
# print "$gene_name\n";
 my $gene_protein = $gene_protein_hash{$gene_id};
    $elegans_protein_id="";
    $ensp="";
    if ($gene_protein){
        $ensp = $species_ensp_hash{$gene_protein};
        $elegans_protein_id = $species_ce_protein_hash{$gene_protein};
    }
    if ($elegans_protein_id){
        $elegans_ortholog_gene = $elegans_name_hash{$elegans_protein_id};
    }
    if ($ensp){
    if (length($ensp) gt 1){
        $ensp =~s/ENSEMBL\://gi;
        $hgnc = $biomart_hgnc_hash{$ensp};
        $description = $biomart_description_hash{$ensp};
        }
    }

    if (($description) and ($description =~/$uncharacterized/)){
         $description = "";
    }
   
    if (($description) and ($description !~/$uncharacterized/) and ($hgnc)){
#         $description =~ s/\[(?:[^\[\]]++|(?0))*+]//g;
         $description =~ s/\[.+?\]//g;
         $description =~ s/^\s+//;
         $description =~ s/\s+$//;
# Rule - remove any reference to kDa 
         $description =~s/\,\s\d*kDa//gi;
         $description =~s/\d*kDa//gi;
         $description =~s/\dkDa//gi;
         $description =~s/\(\d*kDa\)//gi;
         $description =~s/\(\d*kD\)//gi;
         $description =~s/\(\d+.\d+kD\)//gi;
         $description =~s/\(\d+.\d+kDa\)//gi;
         $description =~s/\d+.\d+kDa//gi;
# Rule - apply the same rules regarding family, sufamily, member in description names as C. elegans, create_sentence_homologies.pl.
         if ($description =~/$family4/) {
             $description =~s/$family4//g;
          if ($description =~/(\d+)\s(\S+)$/){
              my $substring1 = $1;
              my $substring2 = $2;
              $description =~s/$substring1\s$substring2/$substring1$substring2/g;
            }
         }
         if (($description =~/$member/) and ($description =~/$class/)){
              $description =~ s/$class//g;
          }
         if ( (($description =~/$family/) and ($description =~/$subfamily/))or
              (($description =~/$family/) and ($description =~/$group/))or
              (($description =~/$family/) and ($description =~/$polypeptide/))or
              (($description =~/$family/) and ($description =~/$member/))or
              (($description =~/$subfamily/) and ($description =~/$group/))or
              (($description =~/$subfamily/) and ($description =~/$polypeptide/))or
              (($description =~/$subfamily/) and ($description =~/$member/))or
              (($description =~/$subfamily2/) and ($description =~/$group/))or
              (($description =~/$subfamily2/) and ($description =~/$polypeptide/))or
              (($description =~/$subfamily2/) and ($description =~/$member/))or
              (($description =~/$subfamily2/) and ($description =~/$family/))or
              (($description =~/$group/) and ($description =~/$polypeptide/))or
              (($description =~/$group/) and ($description =~/$member/))or
              (($description =~/$member/) and ($description =~/$polypeptide/)) or 
              (($description =~/$class/) and ($description =~/$subfamily/))or
              (($description =~/$class/) and ($description =~/$group/))or
              (($description =~/$class/) and ($description =~/$polypeptide/))or
              (($description =~/$class/) and ($description =~/$member/))or
              (($description =~/$class/) and ($description =~/$family/)) ) {

      if ($description !~/$family/){
              $description =~s/$subfamily\s?/ /g;
              $description =~s/$subfamily2\s?/ /g;
              $description =~s/$subfamily3\s?//g;
       } else{
              $description =~s/$subfamily//g;
              $description =~s/$subfamily2//g;
       }
      
      $description =~s/$family\s?/ /g;
      $description =~s/$family2\s?/ /g;
      $description =~s/$subfamily\s?/ /g;
      $description =~s/$subfamily2\s?/ /g;
      $description =~s/$group\s?//g;
      $description =~s/$member\s?//g;
      $description =~s/$polypeptide\s?//g;
      $description =~s/$family3//g;
         if ($description =~/\)(\d+)$/){
              my $substring1 = $1;
              $description =~s/\)$substring1/\) $substring1/g;
          }
        }
     }

       my $mod = 0;
       my $target="";
       if (($description) and ($ensp)){
           my @target_matches=($description =~m/\((.+?)\)/gx);
         if (@target_matches gt 0){
          foreach my $element (@target_matches){
           chomp($element);
           $element =~ s/^\s+//;
           $element =~ s/\s+$//;
           if ( first {$_ eq $element} @mods ){
            $target = $element;
           } 
           if ($element=~/$elegans/){
            $target = $elegans;           
          } 
        }
       }
      if (($target) and (@target_matches) ) {
       if ( first {$_ eq $target} @mods ){
          $mod =1;
# remove the parentheses and what is inside
          $description =~s/\s*\([^)]*\)//g;
# reorder based on rule.
          $description =~s/$homolog//g;
          $description =~s/$homolog2//g;
          my $temp_string = $description;
          $temp_string =~s/\s\,/\,/g;
          $description = $target . $human . $temp_string;
          $description =~ s/^\s+//;
          $description =~ s/\s+$//;
         } 
         if ($target =~/$elegans/){
           $mod =2;
         } 
        }

       }

      if (($elegans_ortholog_gene) and ($description) and ($ensp)){
      $sentence  = $gene_name;
      if (($mod==0) or ($mod==2)){
          $sentence .= $elegans_string . $elegans_ortholog_gene . $and_human_string . $description;
      } elsif ($mod ==1) {
          $sentence .= $elegans_string1 . $elegans_ortholog_gene . $and_of . $description;
      }
      $sentence =~ s/,+$//;
      $sentence .= " $hgnc\;\n";
      write_file($output_file, {append => 1 }, $sentence);
      write_file($output_file, {append => 1 }, "\n");
      write_file($out, $sentence);
     } elsif (($description) and ($ensp) and not ($elegans_ortholog_gene)) {
      $sentence  = $gene_name;
      if (($mod==0) or ($mod==2)){
          $sentence .= $only_human_string . $description;
      } elsif ($mod ==1) {
          $sentence .= $only_human_string1 . $description;
      }
      $sentence =~ s/,+$//;
      $sentence .= " $hgnc\;\n";
      write_file($output_file, {append => 1 }, $sentence);
      write_file($output_file, {append => 1 }, "\n");
      write_file($out, $sentence);
     } elsif (($elegans_ortholog_gene) and not ($description) and not ($ensp)){
      $sentence  = $gene_name;
      $sentence .= $elegans_string . $elegans_ortholog_gene;
      $sentence .= "\;\n";
      write_file($output_file, {append => 1 }, $sentence);
      write_file($output_file, {append => 1 }, "\n");
      write_file($out, $sentence);
    } 
# else {
#      print "$gene_name not included\.\n";
#    }
   
}
exit 0;

sub get_gene_array{
 my $xrefs=shift;
 my @lines = read_file($xrefs);
 my @array=();
 foreach my $line (@lines){
    my ($gene_id0, $wb_gene_id, $gene_name, $gene_id1, $gene_protein_id, $id0, $id1) = split(/\t/, $line);
        $wb_gene_id =~ s/^\s+//;
        $wb_gene_id =~ s/\s+$//;
    push(@array, $wb_gene_id);
 }
 return \@array;
}

sub get_gene_name_hash{
 my $species = shift;
 my $xrefs   = shift;
 my %gene_name_hash=();
 my @lines = read_file($xrefs);
 foreach my $line (@lines){
  my ($gene_id0, $wb_gene_id, $gene_name, $gene_id1, $gene_protein_id, $id0, $id1) = split(/\t/, $line);
        $wb_gene_id =~ s/^\s+//;
        $wb_gene_id =~ s/\s+$//;
        $gene_name =~ s/^\s+//;
        $gene_name =~ s/\s+$//;
        $gene_id0 =~ s/^\s+//;
        $gene_id0 =~ s/\s+$//;
  my $key = $wb_gene_id;
  my $value = $gene_name;
  if ($gene_name eq "\."){
      $value = $gene_id0;
  }
      $gene_name_hash{$key} = $value;
 }
      return \%gene_name_hash;
}
sub get_inverse_gene_protein_hash{
 my $species = shift;
 my $xrefs   = shift;
 my %inverse_gene_protein_hash=();
 my @lines = read_file($xrefs);
 foreach my $line (@lines){
  my ($gene_id0, $wb_gene_id, $gene_name, $gene_id1, $gene_protein_id, $id0, $id1, $id2) = split(/\t/, $line);
        $gene_protein_id =~ s/^\s+//;
        $gene_protein_id =~ s/\s+$//;
        $gene_name =~ s/^\s+//;
        $gene_name =~ s/\s+$//;
        $gene_id0 =~ s/^\s+//;
        $gene_id0 =~ s/\s+$//;

        my $value = $gene_name;
        if (($value eq "\.") or ($value eq "") or not ($value)){
            $value = $gene_id0;
        }
           $value =~ s/^\s+//;
           $value =~ s/\s+$//;
      if ($gene_protein_id){
          $inverse_gene_protein_hash{$gene_protein_id} = $value;
      }
 }
 return \%inverse_gene_protein_hash;
}
sub get_species_ce_protein_blastp_hash{
   my $species = shift;
   my $file  = shift;
   my %hash=();
   my $csv = Text::CSV->new ({
     binary    => 1,
     auto_diag => 1,
     sep_char  => ','    # not really needed as this is the default
   });
   open(my $data, '<:encoding(utf8)', $file) or die "Could not open '$file' $!\n";
   while (my $fields = $csv->getline( $data )) {
   
        my $protein_id = $fields->[0]; 
        my $ce_protein_id = $fields->[19];
        $protein_id =~ s/^\s+//;
        $protein_id =~ s/\s+$//;
#        print "$ce_protein_id\n";
        if ($species =~/pacificus/){
            $ce_protein_id = $fields->[17];
        }
        $ce_protein_id =~ s/WP\://gi;
        $ce_protein_id =~ s/WP\://gi;
        $ce_protein_id =~ s/^\s+//;
        $ce_protein_id =~ s/\s+$//;
        $hash{$protein_id} = $ce_protein_id;
    }
  close($data); 
  return \%hash; 
}
sub get_species_elegans_protein_hash{
   my $species = shift;
   my $file  = shift;
   my %species_elegans_protein_hash=();
   
   my $csv = Text::CSV->new ({
     binary    => 1,
     auto_diag => 1,
     sep_char  => ','    # not really needed as this is the default
   });
   open(my $data, '<:encoding(utf8)', $file) or die "Could not open '$file' $!\n";
   while (my $fields = $csv->getline( $data )) {
        my $protein_id = $fields->[0]; 
        my $elegans_id = $fields->[19];
        $protein_id =~ s/^\s+//;
        $protein_id =~ s/\s+$//;
        $elegans_id =~ s/^\s+//;
        $elegans_id =~ s/\s+$//;
        $species_elegans_protein_hash{$protein_id} = $elegans_id;
    }
  close ($data);
  return \%species_elegans_protein_hash; 
}
sub get_biomart_description_hash{
 my $file = shift;
 my %hash=();
   my $csv = Text::CSV->new ({
     binary    => 1,
     auto_diag => 1,
     sep_char  => ','    # not really needed as this is the default
   });
   open(my $data, '<:encoding(utf8)', $file) or die "Could not open '$file' $!\n";
   while (my $fields = $csv->getline( $data )) {
 
  my $description =  $fields->[2];
  my $ensp = $fields->[8];
  if ($ensp){
      $ensp =~ s/^\s+//;
      $ensp =~ s/\s+$//;
      $ensp =~s/ENSEMBL\://gi;
      $hash{$ensp} = $description;
  }
 }
 close($data);
 return \%hash;
}
sub get_biomart_hgnc_hash{
 my $file = shift;
 my %hash=();

   my $csv = Text::CSV->new ({
     binary    => 1,
     auto_diag => 1,
     sep_char  => ','    # not really needed as this is the default
   });
   open(my $data, '<:encoding(utf8)', $file) or die "Could not open '$file' $!\n";
   while (my $fields = $csv->getline( $data )) {
 
    my $hgnc =  $fields->[4];
    my $ensp = $fields->[8];
       $ensp =~ s/^\s+//;
       $ensp =~ s/\s+$//;
       $hgnc =~ s/^\s+//;
       $hgnc =~ s/\s+$//;
       $ensp =~s/ENSEMBL\://gi;
  if ($ensp){
   $hash{$ensp} = "\(HGNC\:$hgnc\)";
  }
 }
 close($data);
 return \%hash;
}
