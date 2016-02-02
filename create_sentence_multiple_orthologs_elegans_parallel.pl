#!/usr/bin/env perl
use ConciseDescriptions;
use List::MoreUtils qw(uniq);
use List::Util qw/first/;
use File::Slurp;
use Text::CSV;
use Spreadsheet::XLSX;

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
# The list of genes with no concise descriptions is in $infile
my $gene_list_dir = $home . "release/$PRODUCTION_RELEASE/$species/gene_lists/";
my $elegan_gene_list_dir = $home . "release/$PRODUCTION_RELEASE/c_elegans/gene_lists/";
my $db_gene_list  = $elegan_gene_list_dir . "wb_gene_list.txt";

my $orthology = $home . "release/$PRODUCTION_RELEASE/$species/orthology/";
my $old_orthology = $home . "release/$RELEASE/$species/orthology/";
my $elegans_sapiens_orthology = $orthology . "input_files/orthologs.Caenorhabditis_elegans.H_Sapiens.txt";
my $orthology_elegans = $orthology . "input_files/c_elegans.orthologs.txt";
my $output_file = $orthology . "output_files/sentences_for_orthology.txt";
#
if (-e $orthology){
   print "$orthology exists\n"; 
} else {
   mkdir $orthology or die "could not create $orthology";
}
my $output_path = $orthology;
my $individual_path = $output_path . "output_files/individual_gene_sentences/";
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

my $biomart = $orthology . "input_files/HumanIDs_mart_export.txt";
my $reading_biomart_file_for_error = read_file($biomart);
if ($reading_biomart_file_for_error =~ /ERROR/){
    $biomart =  $old_orthology . "input_files/HumanIDs_mart_export.txt";
}
my ($gene_ensg_hash_ref, $gene_hgnc_hash_ref) = get_sapiens_hash($elegans_sapiens_orthology);
my %gene_ensg_hash= %{$gene_ensg_hash_ref};
my %gene_hgnc_hash= %{$gene_hgnc_hash_ref};

my $gene_name_hash_ref = get_gene_name_hash($orthology_elegans, $db_gene_list);
my %gene_name_hash = %$gene_name_hash_ref;

my ($biomart_description_hash_ref, $biomart_hgnc_hash_ref) = get_biomart_description_hash($biomart);
my %biomart_description_hash = %$biomart_description_hash_ref;
my %biomart_hgnc_hash = %$biomart_hgnc_hash_ref;
my $human_gene_families_tsv =  $orthology . "input_files/human_gene_families.txt";
my ($hgnc_family_description_hash_ref,$hgnc_family_hash_ref)=get_human_gene_family_hashes($human_gene_families_tsv);
my $hgnc_gene_family_tag_hash_ref=get_human_gene_family_tag_hash($human_gene_families_tsv);
my %hgnc_family_description_hash = %{$hgnc_family_description_hash_ref};
my %hgnc_family_hash = %{$hgnc_family_hash_ref};
my %hgnc_gene_family_tag_hash = %{$hgnc_gene_family_tag_hash_ref}; 
foreach my $key (keys %hgnc_gene_family_tag_hash) {
    my $value = $hgnc_gene_family_tag_hash{$key};
    print "tag hash $key = \t$value\n";
}
my $and = "and";
my $comma = "\,";
my $is_a = " is an ortholog of human ";
my $human_comma = "human \,";
my $ihuman = "human";
my $doublespace = "  ";
my $space = " ";
my $spacesemicolon = " \;";
my $semicolon ="\;";

# Define the various model organisms used
my @mods=("S\. cerevisiae","yeast","S\. pombe","mouse","Drosophila","Chlamydomonas","zebrafish","Rat","rat","chicken","Strongylocentrotus purpuratus","Arabidopsis","Dictyostelium","bacterial","E\. coli","H\. influenzae","Dictyostelium");
 my $elegans = "C\. elegans";
 my $human = " and human ";
 
my $gene_array_ref = get_gene_array($orthology_elegans);
my @gene_array = @{$gene_array_ref};

my $dead_gene_list = $gene_list_dir . "sort.dead_genes.txt";
my @dead_genes = read_file($dead_gene_list);

my @live_gene_array = ();
foreach my $test (@gene_array){
 my $keep =0;
 foreach my $dead (@dead_genes){
   if ($dead =~/$test/){
    $keep = 1;
   }
 }
 if ($keep ==0){
   push(@live_gene_array, $test);
 }
}
my @sorted_live_gene_array = sort(@live_gene_array);
my $curated_gene_list = $gene_list_dir . "sort.curated_genes.txt";
my @curated_genes = read_file($curated_gene_list);
my @uncurated_genes_array = ();
foreach my $test (@sorted_live_gene_array){
 my $keep =0;
# foreach my $curated (@curated_genes){
#   if ($curated =~/$test/){
#    $keep = 1;
#   }
# }
 if ($keep ==0){
   push(@uncurated_genes_array, $test);
 }
}
my @sorted_uncurated_genes_array = sort(@uncurated_genes_array); 

foreach my $gene_id (@sorted_uncurated_genes_array){
  my $sentence  = "";
  chomp($gene_id);
       $gene_id =~ s/^\s+//;
       $gene_id =~ s/\s+$//;
  my $out = $individual_path . $gene_id;
  my $gene_name = $gene_name_hash{$gene_id};
  my $ensg_list = "";
  my @ensg_array = ();
  if ($gene_ensg_hash{$gene_id}){
      $ensg_list = $gene_ensg_hash{$gene_id};
   if ($ensg_list =~/AND/){
      @ensg_array = split(/ AND /, $ensg_list);
    } else {
     $ensg_array[0] = $ensg_list;
   }
  }
  my $ensg_phrase = "";
  my $count = 0;
  my $ensg_size = @ensg_array;
  if ($ensg_size <= 3){
  foreach my $ensg (@ensg_array){
       $ensg =~ s/^\s+//;
       $ensg =~ s/\s+$//;
    my $gene_hgnc = $gene_hgnc_hash{$ensg};
    my $description = $biomart_description_hash{uc $gene_hgnc};
    if ($gene_hgnc){
    if ($count == 0){
     if ($description){
      $ensg_phrase = $gene_hgnc . " \(" . $description . "\)";  
     } else {
      $ensg_phrase = $gene_hgnc;  
     }
    } elsif (($count > 0) and ($count < ($ensg_size -1))){
     if ($description){
      $ensg_phrase .= "\, " . $gene_hgnc . " \(" . $description . "\)";  
     } else {
      $ensg_phrase .= "\, " . $gene_hgnc;  
     }
    } elsif (($count > 0) and ($count == ($ensg_size -1))){
     if ($description){
      $ensg_phrase .= " and " . $gene_hgnc . " \(" . $description . "\)";
     } else {
      $ensg_phrase .= " and " . $gene_hgnc;
     }  
    } 
   }
    $count++;
  } # end foreach my $ensg (@ensg_array)
     $sentence = $gene_name . $is_a . $ensg_phrase;
  } else {
     my @human_gene_array = ();
     foreach my $ensg (@ensg_array) {
        my $gene_hgnc = $gene_hgnc_hash{$ensg};
        push(@human_gene_array, $gene_hgnc);
     }
#  while( my( $key, $value ) = each %hgnc_gene_family_tag_hash ){
#    print "before reduce orthologs hgnc gene family tag hash $key\t$value\n";
#}

     $ensg_phrase = reduce_orthologs(\@human_gene_array,
                                     \%hgnc_family_description_hash,
                                     \%hgnc_family_hash,
                                     \%hgnc_gene_family_tag_hash,
                                     \%biomart_description_hash);
     $sentence = $gene_name . $ensg_phrase;

  }

  if ($ensg_phrase){
#      $sentence = $gene_name . $is_a . $ensg_phrase;
#      $sentence =~ s/$human_comma/$ihuman/g;
      $sentence =~ s/^\s+//g;
      $sentence =~ s/\s+$//g;
      $sentence =~ s/$doublespace//g;

      $sentence .= "\;\n";
      write_file($output_file, {append => 1 }, $gene_id);
      write_file($output_file, {append => 1 }, "\n");
      write_file($output_file, {append => 1 }, $sentence);
      write_file($output_file, {append => 1 }, "\n\n\n");
      write_file($out, $sentence);
 }
 
} # end foreach my $gene_id (@gene_array)
exit 0;
sub get_gene_array{
 my $orthology = shift;
 my @gene_array=();
 my @lines = read_file($orthology);
 foreach my $line (@lines){
  next if ($line =~ /\#/);
  next if ($line =~ /\=/);
  if ($line =~/WBGene/){
     my ($WBGene, $Gene) = split(/\t/, $line);         
        $WBGene =~ s/^\s+//;
        $WBGene =~ s/\s+$//;
        $Gene =~ s/^\s+//;
        $Gene =~ s/\s+$//;
        push(@gene_array, $WBGene);
     }
 }
      my @sorted_gene_array = uniq(sort(@gene_array));
      return \@sorted_gene_array;
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

         if ($multiple_size > 0){
             foreach my $a_line (@multiple_source_array){
               my ($a_wb_gene_id, $a_gene_species_name, $a_ensg, $a_hgnc, $a_source) = split(/\t/, $a_line);
               my ($a_gene_name, $a_species_name) = split(/ /, $a_gene_species_name);
                   $gene_hgnc_hash{$a_ensg}  = $a_hgnc;

               if ($gene_ensg_hash{$a_wb_gene_id}){
                   $gene_ensg_hash{$a_wb_gene_id} .= " AND " .  $a_ensg;
               } else{
                   $gene_ensg_hash{$a_wb_gene_id} = $a_ensg;
               }
            }
          } else {
             foreach my $a_line (@single_source_array){
               my ($a_wb_gene_id, $a_gene_species_name, $a_ensg, $a_hgnc, $a_source) = split(/\t/, $a_line);
               my ($a_gene_name, $a_species_name) = split(/ /, $a_gene_species_name);
                   $gene_hgnc_hash{$a_ensg}  = $a_hgnc;

               if ($gene_ensg_hash{$a_wb_gene_id}){
                   $gene_ensg_hash{$a_wb_gene_id} .= " AND " .  $a_ensg;
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
          if ($multiple < 50){
           push(@multiple_source_array, $line);
          }
       } else {
          my $single = @single_source_array;
          if ($single < 50){
           push(@single_source_array, $line);
          }
       }
  }
 }
 return (\%gene_ensg_hash, \%gene_hgnc_hash);
}
 
sub get_gene_name_hash{
 my $orthology = shift;
 my $db_gene_list = shift;
 my %gene_name_hash=();
 my @lines = read_file($orthology);
 foreach my $line (@lines){
  next if ($line =~ /\#/);
  next if ($line =~ /\=/);
  if ($line =~/WBGene/){
     my ($WBGene, $Gene) = split(/\t/, $line);         
        $WBGene =~ s/^\s+//;
        $WBGene =~ s/\s+$//;
        $Gene =~ s/^\s+//;
        $Gene =~ s/\s+$//;
        $gene_name_hash{$WBGene} = $Gene;
     }
 }
  my @db_lines = read_file($db_gene_list);
 foreach my $db_line (@db_lines){
  chomp($db_line);
  if ($db_line =~/WBGene/){
     my ($WBGene, $Gene) = split(/\t/, $db_line);         
        $WBGene =~ s/^\s+//;
        $WBGene =~ s/\s+$//;
        $Gene =~ s/^\s+//;
        $Gene =~ s/\s+$//;
        $gene_name_hash{$WBGene} = $Gene;
     }
 }
      return \%gene_name_hash;
}
sub get_biomart_description_hash{
 my $file = shift;
 my %hash=();
 my %ahash=();

# define some variables
my $family = "\, family ";
my $family2 = " family member";
my $family3 = " family";
my $family4 = " family\, member";
my $uncharacterized = "uncharacterized";
my $group = "\, group ";
my $member = "\, member ";
my $subfamily = "\, subfamily ";
my $subfamily2 = "\, sub\-family ";
my $subfamily3 = "subfamily";
my $polypeptide ="\, polypeptide ";
my $class = "\, class ";
my $homolog  = "homolog ";
my $homolog2 = "homolog";
my $homolog3 = "homolog\.*";
#
   my $csv = Text::CSV->new ({
     binary    => 1,
     auto_diag => 1,
     sep_char  => ','    # not really needed as this is the default
   });
   open(my $data, '<:encoding(utf8)', $file) or die "Could not open '$file' $!\n";
   while (my $fields = $csv->getline( $data )) {
 
  my $description =  $fields->[2];
  my $ensg = $fields->[0];
  my $hgnc = uc($fields->[4]);
       $ensg =~ s/^\s+//;
       $ensg =~ s/\s+$//;
       $hgnc =~ s/^\s+//;
       $hgnc =~ s/\s+$//;
       
  if ($ensg){
      $ensg =~s/ENSEMBL\://gi;
      if ($description !~/$uncharacterized/i){
         $description =~ s/\[.+?\]//g;
         $description =~ s/^\s+//;
         $description =~ s/\s+$//;
         $description =~s/$homolog//g;
         $description =~s/$homolog2//g;
# Rule - remove any reference to kDa 
         $description =~s/\,\s\d*kDa//gi;
         $description =~s/\d*kDa//gi;
         $description =~s/\dkDa//gi;
         $description =~s/\(\d*kDa\)//gi;
         $description =~s/\(\d*kD\)//gi;
         $description =~s/\(\d+.\d+kD\)//gi;
         $description =~s/\(\d+.\d+kDa\)//gi;
         $description =~s/\d+.\d+kDa//gi;
# Rule - apply the same rules regarding family, subfamily, member in description names 
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
          $description =~s/ \(\)//g;  
          $hgnc        =~ s/^\s+//;
          $hgnc        =~ s/\s+$//;
          $description =~ s/^\s+//;
          $description =~ s/\s+$//;

      if (($hgnc) and ($description)){
          $hash{$hgnc} = $description;
          my $key = quotemeta($description);
          $ahash{$key} = $hgnc;
      }
      }
  }
 }
 close($data);
 return (\%hash, \%ahash);
}
sub get_human_gene_family_tag_hash{
my $filename=shift;
my @tsv_lines = read_file($filename, binmode => ':utf8');
my %hgnc_gene_family_tag_hash=();

foreach my $line (@tsv_lines){
 next if ($line =~/HGNC/);
 chomp($line);
 my @elements = split(/\t/, $line);
 my $hgnc_symbol = $elements[1];
 my $hgnc_name   = $elements[2];
 my $hgnc_gene_family = $elements[9];
 my $hgnc_gene_family_description = $elements[10];

    $hgnc_name   =~ s/^\s+//;
    $hgnc_name   =~ s/\s+$//;
    $hgnc_symbol =~ s/^\s+//;
    $hgnc_symbol =~ s/\s+$//;
    $hgnc_gene_family =~ s/^\s+//;
    $hgnc_gene_family =~ s/\s+$//;
    $hgnc_gene_family_description =~ s/^\s+//;
    $hgnc_gene_family_description =~ s/\s+$//;
    my $location = index ( $hgnc_gene_family_description, "," );
    if ( $location > -1 ){
         $hgnc_gene_family_description =~s/\, family\s\d+//gi;
         $hgnc_gene_family_description =~s/\,family\s\d+//gi;
         $hgnc_gene_family_description =~s/family//gi;
         $hgnc_gene_family_description =~s/superfamily//gi;
         $hgnc_gene_family_description =~ s/^\s+//;
         $hgnc_gene_family_description =~ s/\s+$//;
     my @parts = split(/\,/, $hgnc_gene_family_description);
     if (@parts == 2){
         $parts[0] =~ s/^\s+//;
         $parts[0] =~ s/\s+$//;
         $parts[1] =~ s/^\s+//;
         $parts[1] =~ s/\s+$//;
       $hgnc_gene_family_description = $parts[1] . " " . $parts[0];
     }
    }

 if ($hgnc_gene_family){
  if ($hgnc_gene_family_description){
   my $key = quotemeta($hgnc_gene_family_description); 
   $hgnc_gene_family_tag_hash{$key} = $hgnc_gene_family;
  }
 }
} 

foreach my $key (keys %hgnc_gene_family_tag_hash) {
    my $value = $hgnc_gene_family_tag_hash{$key};
    print "sub $key\t$value";
}

return \%hgnc_gene_family_tag_hash;
}
sub get_human_gene_family_hashes{
my $filename=shift;
my @tsv_lines = read_file($filename, binmode => ':utf8');
my %hgnc_family_hash=();
my %hgnc_family_description_hash=();
my %hgnc_family_tag_hash=();

foreach my $line (@tsv_lines){
 next if ($line =~/HGNC/);
 chomp($line);
 my @elements = split(/\t/, $line);
 my $hgnc_symbol = $elements[1];
 my $hgnc_name   = $elements[2];
 my $hgnc_gene_family = $elements[9];
 my $hgnc_gene_family_description = $elements[10];

    $hgnc_name   =~ s/^\s+//;
    $hgnc_name   =~ s/\s+$//;
    $hgnc_symbol =~ s/^\s+//;
    $hgnc_symbol =~ s/\s+$//;
    $hgnc_gene_family =~ s/^\s+//;
    $hgnc_gene_family =~ s/\s+$//;
    $hgnc_gene_family_description =~ s/^\s+//;
    $hgnc_gene_family_description =~ s/\s+$//;
 $hgnc_family_description_hash{$hgnc_symbol} = $hgnc_gene_family_description;
 if ($hgnc_gene_family){
  if ($hgnc_symbol){
   $hgnc_family_hash{$hgnc_symbol} = $hgnc_gene_family;
  }
 }
} 
return (\%hgnc_family_description_hash, \%hgnc_family_hash);
}

sub reduce_orthologs{
 my $human_gene_array_ref = shift;
 my $hgnc_family_description_hash_ref = shift;
 my $hgnc_family_hash_ref = shift;
 my $hgnc_gene_family_tag_hash_ref = shift;
 my $description_hash_ref = shift;

 my @human_gene_array = @{$human_gene_array_ref};
 my %hgnc_family_description_hash = %{$hgnc_family_description_hash_ref};
 my %hgnc_family_hash = %{$hgnc_family_hash_ref}; 
 my %hgnc_gene_family_tag_hash = %{$hgnc_gene_family_tag_hash_ref};
 my %description_hash = %{$description_hash_ref};

 my $phrase = "";
 my %family_hash=();
 my %family_elements_hash=();
 my @families=();
 my @unclassified=();

  while( my( $key, $value ) = each %hgnc_gene_family_tag_hash ){
    print "hgnc gene family tag hash $key\t$value\n";
}

  foreach my $ortholog (@human_gene_array){
   my $family="";
   print "ortholog is $ortholog\n";
   if ($hgnc_family_description_hash{$ortholog}){
       $family = $hgnc_family_description_hash{$ortholog};
   } else {
       $family = "unclassified";
       my $txt_description=$description_hash{$ortholog};
       print "$ortholog\t$txt_description is not classified\n";
       push(@unclassified, $ortholog);
   }
   push(@families, $family);
   my $key = quotemeta($family);
   if ($family_hash{$key}){
     $family_elements_hash{$key}+= 1; 
     $family_hash{$key} .= "\,$ortholog";
   } else {
     $family_elements_hash{$key} = 1; 
     $family_hash{$key} = $ortholog;
   }   
 }

 @families = uniq(@families);
 @unclassified = uniq(@unclassified);

 my @singlets=();
 my @multiples=();

foreach my $f (@families){
 my $count=0;
 next if ($f=~/unclassified/);
   my $key = quotemeta($f);
   $count = $family_elements_hash{$key};
   if ($count == 1){
    push (@singlets, $f);
   } elsif ($count > 1) {
    push (@multiples, $f);
  }
 }
 @singlets = uniq(@singlets);
 @multiples = uniq(@multiples);

my $u_size=@unclassified;
my $s_size=0;
foreach my $s (@singlets){
  next if ($s =~/unclassified/);
  $s_size++;
}
my $m_size=0;
foreach my $m (@multiples){
  next if ($m =~/unclassified/);
  $m_size++;
}
my $kount=0;
if (($s_size > 0) and ($u_size==0)) {
 $phrase = " is an ortholog of human ";
foreach my $s (@singlets){
      $s =~s/^\s+//;
      $s =~s/\s+$//;
   my $key = quotemeta($s);
   my $gene_list = $family_hash{$key};
   my @genes = split(/\,/, $gene_list);
 foreach my $gene (@genes) {
   my $description = $description_hash{uc $gene};
  if ($kount == 0){
     if ($description){
      $phrase .= $gene . " \(" . $description . "\)";
     } else {
      $phrase .= $gene;
     }
    } elsif (($kount > 0) and ($kount < ($s_size-1))){
     if ($description){
      $phrase .= "\, " . $gene . " \(" . $description . "\)";
     } else {
      $phrase .= "\, $gene";
     }  
    } elsif (($kount > 0) and ($kount == ($s_size-1)) and ($s_size > 1)){
     if ($description){
      $phrase .= " and " . $gene . " \(" . $description . "\)";
     } else {
      $phrase .= " and $gene";
     }    
    }
   $kount++;
  }
 }
}
$kount=0;
if (($s_size > 0) and ($u_size > 0)) {
 $phrase = " is an ortholog of human ";
foreach my $s (@singlets){
      $s =~s/^\s+//;
      $s =~s/\s+$//;
   my $key = quotemeta($s);
   my $gene_list = $family_hash{$key};
   my @genes = split(/\,/, $gene_list);
 foreach my $gene (@genes) {
   my $description = $description_hash{uc $gene};
  if ($kount == 0){
     if ($description){
      $phrase .= $gene . " \(" . $description . "\)";
     } else {
      $phrase .= $gene;
     }
    } elsif (($kount > 0) and ($kount < $s_size)){
     if ($description){
      $phrase .= "\, " . $gene . " \(" . $description . "\)";
     } else {
      $phrase .= "\, $gene";
     }  
    } 
   $kount++;
   }
  }
   foreach my $uc (@unclassified){
   my $description = $description_hash{uc $uc};
  if ($kount == 0){
     if ($description){
      $phrase .= $uc . " \(" . $description . "\)";
     } else {
      $phrase .= $uc;
     }
    } elsif (($kount > 0) and ($kount < ($u_size-1))){
     if ($description){
      $phrase .= "\, " . $uc . " \(" . $description . "\)";
     } else {
      $phrase .= "\, $uc";
     }  
    } elsif (($u_size > 1) and ($kount == ($u_size-1))){
     if ($description){
      $phrase .= " and " . $uc . " \(" . $description . "\)";
     } else {
      $phrase .= " and $uc";
     } 
    }
   $kount++;
  }
 }
$kount=0;
if (($u_size > 0) and ($s_size==0)) {
 $phrase = " is an ortholog of human ";
   foreach my $uc (@unclassified){
   my $description = $description_hash{uc $uc};
  if ($kount == 0){
     if ($description){
      $phrase .= $uc . " \(" . $description . "\)";
     } else {
      $phrase .= $uc;
     }
    } elsif (($kount > 0) and ($kount < ($u_size-1))){
     if ($description){
      $phrase .= "\, " . $uc . " \(" . $description . "\)";
     } else {
      $phrase .= "\, $uc";
     }  
    } elsif (($u_size > 1) and ($kount == ($u_size-1))){
     if ($description){
      $phrase .= " and " . $uc . " \(" . $description . "\)";
     } else {
      $phrase .= " and $uc";
     } 
    }
   $kount++;
  }
 }
$kount = 0;
if ($m_size > 0){
 if ($phrase){
  $phrase .= " and members of the ";
 } else {
  $phrase = " is an ortholog of members of the human ";
 }
foreach my $m (@multiples){
      $m =~s/^\s+//;
      $m =~s/\s+$//;
   my $key = quotemeta($m);
   my $gene_list = $family_hash{$key};
   my $hgnc = $hgnc_gene_family_tag_hash{$key};
      print "hgnc is $hgnc\n";
   my @genes = split(/\,/, $gene_list);
   my $first_gene = $genes[0];
   print "first gene is $first_gene\n";
      $m =~s/superfamily//gi;
      $m =~s/family//gi;
      $m =~s/^\s+//;
      $m =~s/\s+$//;
  if ($kount == 0){
     if ($hgnc){
      $phrase .= $hgnc . " \(" . $m . "\) family including " . $first_gene;
     } else {
      $phrase .= $m . " family including " . $first_gene;
     }
    } elsif (($kount > 0) and ($kount < ($m_size-1))){
     if ($hgnc){
      $phrase .= "\, " . $hgnc . " \(" . $m . "\) family including " . $first_gene;
     } else {
      $phrase .= "\, $m family including " . $first_gene;
     }  
    } elsif (($kount > 0) and ($kount == ($m_size-1)) and ($m_size > 1)){
     if ($hgnc){
      $phrase .= " and " . $hgnc . " \(" . $m . "\) family including " . $first_gene;
     } else {
      $phrase .= " and $m family including " . $first_gene;
     }    
    }
   $kount++; 
  }
 }  
 return $phrase;
}
