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
my $name = $ARGV[2];
my $prefix= $ARGV[3];

chomp($species);
chomp($project);
chomp($name);
chomp($prefix);
       
$project =~ s/^\s+//;
$project =~ s/\s+$//;
$species =~ s/^\s+//;
$species =~ s/\s+$//;
$name    =~ s/^\s+//;
$name    =~ s/\s+$//;
$prefix  =~ s/^\s+//;
$prefix  =~ s/\s+$//;

my $RELEASE = ConciseDescriptions::get_release();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
# The list of genes with no concise descriptions is in $infile
my $elegans_gene_list_dir = $home . "release/$PRODUCTION_RELEASE/c_elegans/gene_lists/";
my $db_gene_list  = $elegans_gene_list_dir . "wb_gene_list.txt";

my $curated_gene_list = $elegans_gene_list_dir . "sort.curated_genes.txt";
my @curated_genes = read_file($curated_gene_list);

my $orthology = $home . "release/$PRODUCTION_RELEASE/$species/orthology/";
my $orthology_elegans_dir = $home . "release/$PRODUCTION_RELEASE/c_elegans/orthology/";

my $sapiens_orthology = $orthology . "input_files/orthologs.$name.H_Sapiens.txt";
my $orthology_elegans = $orthology_elegans_dir . "input_files/c_elegans.orthologs.txt";
my $orthology_elegans_species = $orthology . "input_files/orthologs.$species.elegans.txt";
my $orthology_species = $orthology . "input_files/$species.orthologs.txt";
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

my ($gene_elegans_id_hash_ref, $gene_elegans_hash_ref) = get_id_name_hash($orthology_elegans_species, $db_gene_list);
my %gene_elegans_hash= %{$gene_elegans_hash_ref};
my %gene_elegans_id_hash= %{$gene_elegans_id_hash_ref};

my $gene_name_hash_ref = get_gene_name_hash($orthology_species, $db_gene_list);
my %gene_name_hash = %$gene_name_hash_ref;

my $gene_array_ref = get_gene_array($orthology_species);
my @gene_array = @{$gene_array_ref};

my $and = "and";
my $AND = "AND";
my $comma = "\,";
my $is_a = " is an ortholog of human ";
my $is_a_ce = " is an ortholog of C\. elegans ";
my $human = " and an ortholog of human ";

my $doublespace = "  ";
my $space = " ";
my $spacesemicolon = " \;";
my $semicolon = "\;";

# Define the various model organisms used
my @mods=("S\. cerevisiae","yeast","S\. pombe","mouse","Drosophila","Chlamydomonas","zebrafish","Rat","rat","chicken","Strongylocentrotus purpuratus","Arabidopsis","Dictyostelium","bacterial","E\. coli","H\. influenzae","Dictyostelium");
 my $elegans = "C\. elegans";

my $dead_gene_list = $elegans_gene_list_dir . "sort.dead_genes.txt";
my @dead_genes = read_file($dead_gene_list);

my @uncurated_genes_array = ();
foreach my $test (@gene_array){
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
my @uncurated_live_genes_array = ();
foreach my $test (@sorted_uncurated_genes_array){
 my $keep =0;
 foreach my $dead (@dead_genes){
   if ($dead =~/$test/){
    $keep = 1;
   }
 }
 if ($keep ==0){
   push(@uncurated_live_genes_array, $test);
 }
}
my @sorted_uncurated_live_genes_array = uniq(sort(@uncurated_live_genes_array)); 
 
foreach my $gene_id (@sorted_uncurated_live_genes_array){
  my $sentence  = "";
  chomp($gene_id);
       $gene_id =~ s/^\s+//;
       $gene_id =~ s/\s+$//;
  my $out = $individual_path . $gene_id;
  my $gene_name = $gene_name_hash{$gene_id};
#
# Fill array, @elegans_array, with list of WBGeneID orthologs.
#
  my $elegans_list = "";
  my @elegans_array = ();
  if ($gene_elegans_id_hash{$gene_id}){
      $elegans_list = $gene_elegans_id_hash{$gene_id};
   if ($elegans_list =~/$AND/){
      @elegans_array = split(/$AND/, $elegans_list);
    } else {
     $elegans_array[0] = $elegans_list;
   }
  }
#
# Now find the GO process terms associated with each element in the @elegans_array.
#
  my $elegans_phrase = "";
  my $count_elegans = 0;
  my $elegans_size = @elegans_array;
  foreach my $elegans (@elegans_array){
       $elegans =~ s/^\s+//;
       $elegans =~ s/\s+$//;
    my $gene_elegans = $gene_elegans_hash{$elegans};
    if ($gene_elegans){
    if ($count_elegans == 0){
     $elegans_phrase = $gene_elegans;  
    } elsif (($count_elegans > 0) and ($count_elegans < ($elegans_size -1))){
     $elegans_phrase .= "\, " . $gene_elegans;  
    } elsif (($count_elegans > 0) and ($count_elegans == ($elegans_size -1))){
     $elegans_phrase .= " and " . $gene_elegans;  
    }
   }
    $count_elegans++;
  } # end foreach my $elegans

if ($elegans_phrase) {
      $sentence  = $gene_name . $is_a_ce . $elegans_phrase;
      $sentence .= "\;\n";
      write_file($output_file, {append => 1 }, $sentence);
      write_file($output_file, {append => 1 }, "\n\n");
      write_file($out, $sentence);
} 
 
} # end foreach my $gene_id (@gene_array)
exit 0;

sub get_id_name_hash{
 my $orthology=shift;
 my $db_gene_list  = shift;

 my @lines = read_file($orthology);
 my $old_gene="";
 my %gene_name_hash=();
 my %gene_id_hash=();
 my @multiple_source_array=();
 my @single_source_array=();
 my $oldgene="";
 my $newgene; 
 foreach my $line (@lines){ 
    chomp($line);   
    if ($line =~/WBGene/){
     my ($wb_gene_id, $gene_name, $species_name, $id, $name, $source) = split(/\t/, $line); 
        $wb_gene_id =~ s/^\s+//;
        $wb_gene_id =~ s/\s+$//;
        $newgene = $wb_gene_id;
     if ($newgene ne $oldgene){
         my $single_size = @single_source_array;
         my $multiple_size = @multiple_source_array;

         if ($multiple_size > 0){
             foreach my $a_line (@multiple_source_array){
               my ($a_wb_gene_id, $a_gene_species_name, $a_id, $a_name, $a_source) = split(/\t/, $a_line);
               my ($a_gene_name, $a_species_name) = split(/ /, $a_gene_species_name);
                       $gene_name_hash{$a_id}  = $a_name;
 
               if ($gene_id_hash{$a_wb_gene_id}){
                   $gene_id_hash{$a_wb_gene_id} .= "AND" .  $a_id;
               } else{
                   $gene_id_hash{$a_wb_gene_id} = $a_id;
               }
            }
          } elsif ($single_size > 0) {
             foreach my $a_line (@single_source_array){
               my ($a_wb_gene_id, $a_gene_species_name, $a_id, $a_name, $a_source) = split(/\t/, $a_line);
               my ($a_gene_name, $a_species_name) = split(/ /, $a_gene_species_name);
#                   $gene_name_hash{$a_id}  = $a_name;
#                   if not ($gene_name_hash{$a_id}){
                       $gene_name_hash{$a_id}  = $a_name;
#                    } 

               if ($gene_id_hash{$a_wb_gene_id}){
                   $gene_id_hash{$a_wb_gene_id} .= "AND" .  $a_id;
               } else{
                   $gene_id_hash{$a_wb_gene_id} = $a_id;
               }
          }
         }
         @multiple_source_array=();
         @single_source_array=();
         $oldgene=$newgene;
     }
         
         if ($line =~/\;/){
          my $multiple = @multiple_source_array;
           push(@multiple_source_array, $line);
       } elsif ($line) {
          my $single = @single_source_array;
           push(@single_source_array, $line);
       }
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

 return \%gene_id_hash, \%gene_name_hash;
}
sub get_gene_array{
 my $orthology = shift;
 my @gene_array=();
 my @lines = read_file($orthology);
 foreach my $line (@lines){
  next if ($line =~ /\#/);
  next if ($line =~ /\=/);
  chomp($line);
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
sub get_gene_name_hash{
 my $orthology = shift;
 my $db_gene_list = shift;
 my %gene_name_hash=();
 my @lines = read_file($orthology);
 foreach my $line (@lines){
  next if ($line =~ /\#/);
  next if ($line =~ /\=/);
  chomp($line);
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
