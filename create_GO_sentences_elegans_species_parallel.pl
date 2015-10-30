#!/usr/bin/env perl
use ConciseDescriptions;
use LWP::Simple;
use List::MoreUtils qw(uniq);
use List::Util qw/first/;
use File::Slurp;
use Text::CSV;
use Switch;
use warnings;
use strict;

my $AND = "AND";

my $species_term = $ARGV[0];
my ($species, $project, $species_name, $species_prefix) = split(/$AND/, $species_term);

chomp($species);
chomp($project);
chomp($species_name);
chomp($species_prefix);

$project =~ s/^\s+//;
$project =~ s/\s+$//;
$species =~ s/^\s+//;
$species =~ s/\s+$//;
$species_prefix =~ s/^\s+//;
$species_prefix =~ s/\s+$//;
$species_name =~ s/^\s+//;
$species_name =~ s/\s+$//;
$species_name =~ s/\_/ /g;

if ($species !~/elegans/){
#
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
my $html = ConciseDescriptions::get_html_dir();
my $home_elegans = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/c_elegans/";
my $home = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/$species/";
my $gene_list_dir = $home_elegans . "gene_lists/";

my $gene_class_file = $gene_list_dir . "acedb_gene_class.txt";
my $gene_class_hash_ref = ConciseDescriptions::get_gene_class_hash($gene_class_file);
my %gene_class_hash = %{$gene_class_hash_ref};
my $db_gene_list  = $gene_list_dir . "wb_gene_list.txt";
if (-e $db_gene_list){
 print "dp gene list is $db_gene_list\n";
} else {
 print "$db_gene_list does not exist\.\n";
} 
my $orthology = $home . "orthology/";
my $orthology_elegans_species = $orthology . "input_files/orthologs.$species.elegans.txt";
my $orthology_species = $orthology . "input_files/$species.orthologs.txt";
# Define GO term file
my $ontology = "GO";
my $go_file = $gene_list_dir . "go_terms.txt";
my $go_altid =  $gene_list_dir . "goid_altid.txt";
my $input_ace_file = $gene_list_dir . "go_terms.ace";
#
# Create a file containing the GO terms and their respective parents
# build parent/child hash for granularity
#
my ($parents_ref, $children_ref) =  ConciseDescriptions::get_ontology_parents_children($input_ace_file, $ontology);
my %parents = %$parents_ref;
my %children = %$children_ref;
# Create GO hash so that terms can be referenced by GO ID.
my $gene_ontology_ref = ConciseDescriptions::get_ontology_hash($ontology, $go_file, $go_altid);
my %gene_ontology = %$gene_ontology_ref;
my $go_elegans_dir = $home_elegans . "gene_ontology/input_files/";
my $gene_association_elegans_file = $go_elegans_dir . "gene_association.wb";
my $go_species_dir = $home . "gene_ontology/input_files/";
my $gene_association_species_file = $go_species_dir . "gene_association.wb";
my $elegans_orthology = $home_elegans . "orthology/";
my $go_path = $home . "orthology/output_files/";
my $output_file = $go_path . "sentences_from_GO_elegans.txt";
#
my $wbgene_elegans_id_hash_ref = get_wbgene_id_elegans_hash($gene_association_elegans_file, $db_gene_list);
my %wbgene_elegans_id_hash = %$wbgene_elegans_id_hash_ref;
#
my ($gene_elegans_id_hash_ref, $gene_elegans_hash_ref) = get_id_name_hash($orthology_elegans_species, $db_gene_list);
my %gene_elegans_hash= %{$gene_elegans_hash_ref};
my %gene_elegans_id_hash= %{$gene_elegans_id_hash_ref};
#
my $in_c_elegans = " in C\. elegans\, ";
my $doublecomma = "\,\,";
my $the_the = "the the";
my $the = "the";
my $is_is = "is is";
my $is = "is";
my $is_a_is_a = "is a is a";
my $is_a = "is a";
my $is_a_ce = " is an ortholog of C\. elegans ";
my $blank_comma = "\, \,";
my $comma = "\, ";
my $space_comma = " \,";
my $just_comma = "\,";
my $nadp =" NAD or NADP as acceptor\,";
my $component_watch_string = "integral component of";
my $component_watch_string_the = "integral component of the";
my $is_localized_cellular="is localized to the intracellular";
my $is_cellular="is intracellular";
my $structural="structural constituent";
my $is_structural="is a structural constituent";
#
my $synaptic_1_0 = "synaptic transmission\, GABAergic";
my $synaptic_1_1 = "GABAergic synaptic transmission";

my $synaptic_2_0 = "synaptic transmission\, cholinergic";
my $synaptic_2_1 = "cholinergic synaptic transmission";

my $synaptic_3_0 = "synaptic transmission\, dopaminergic";
my $synaptic_3_1 = "dopaminergic synaptic transmission";

my $synaptic_4_0 = "synaptic transmission\, glutamatergic";
my $synaptic_4_1 = "glutamatergic synaptic transmission";

my $synaptic_5_0 = "synaptic transmission\, glycinergic";
my $synaptic_5_1 = "glycinergic synaptic transmission";

my $molting =   "the molting cycle";
my $molting_1 = "molting cycle\, chitin\-based cuticle";
my $molting_2 = "molting cycle\, collagen and cuticulin\-based cuticle";
my $molting_3 = "molting cycle\, protein\-based cuticle";
my $molting_0 = "molting cycle";

my $embryo   = "embryo development";
my $embryo_1 = "embryo development ending in birth or egg hatching";

my $growth_0 = "multicellular organism growth";
my $growth_1 = "growth";

my $the_cell_0 = "localized to the cell\;";
my $the_cell_1 = "expressed widely\;";
#
if (-e $go_path){
   my $go_path_flag = 1;
#   print "$go_path exists\n"; 
} else {
   mkdir $go_path or die "could not create $go_path";
}
my $individual_path = $go_path . "individual_gene_sentences_from_GO_elegans/";
if (-e $individual_path){
   my $individual_path_flag = 1;
#   print "$individual_path exists\n"; 
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

my $gene_array_ref = get_gene_array($orthology_species);
my @gene_array = @{$gene_array_ref};
my $dead_gene_list = $gene_list_dir . "sort.dead_genes.txt";
my @dead_genes = read_file($dead_gene_list);
my $curated_gene_list = $gene_list_dir . "sort.curated_genes.txt";
my @curated_genes = read_file($curated_gene_list);
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
my @sorted_uncurated_live_genes_array = sort(@uncurated_live_genes_array); 

my $gene_name_hash_ref = get_gene_name_hash($orthology_species, $db_gene_list);
my %gene_name_hash = %{$gene_name_hash_ref};
my $go = "P";
my $elegans_process_hash_ref = get_gene_go_hash($gene_association_elegans_file, $go, 1);
my %elegans_process_hash = %{$elegans_process_hash_ref};

foreach my $gene_id (@sorted_uncurated_live_genes_array){
 chomp($gene_id);
       $gene_id =~ s/^\s+//;
       $gene_id =~ s/\s+$//;

 my $gene_name = $gene_name_hash{$gene_id};

 my @gene_processes = ();
 my @GO_processes  = ();
 my $process_string = "";
#
# Obtain the list of orthologs for each gene
#
  my $elegans_list = "";
  my @elegans_array = ();
  my @initial_elegans_array = ();
  if ($gene_elegans_id_hash{$gene_id}){
      $elegans_list = $gene_elegans_id_hash{$gene_id};
   if ($elegans_list =~/$AND/){
      @initial_elegans_array = uniq(split(/$AND/, $elegans_list));
    } else {
     $initial_elegans_array[0] = $elegans_list;
   }
  }
  foreach my $i (@initial_elegans_array){
   chomp($i);
       $i =~ s/^\s+//;
       $i =~ s/\s+$//;
   if ($gene_class_hash{$i}){
    push(@elegans_array, $i);
   } 
  }
  @elegans_array = uniq(sort(@elegans_array));
  my $number_of_orthologs = @elegans_array;  
  my $ortholog_list = "";
  my $kount = 0;
  if ($number_of_orthologs > 0){
    foreach my $element_id (@elegans_array){
     $element_id =~ s/^\s+//;
     $element_id =~ s/\s+$//;
     $kount++;
     my $element = $wbgene_elegans_id_hash{$element_id};
     if ($element){
      if ($kount < ($number_of_orthologs-1)){
       $ortholog_list .= $element . "\, "; 
      } elsif ($kount == ($number_of_orthologs-1)) {
       $ortholog_list .= $element; 
      } elsif (($kount == $number_of_orthologs) and ($number_of_orthologs > 1)) {
       $ortholog_list .= " and " . $element; 
      } elsif ($number_of_orthologs == 1) {
       $ortholog_list = $element; 
      }
     }
    }
  }
  my $ortholog_sentence = "";
  if ($ortholog_list){
      $ortholog_sentence = $in_c_elegans . $ortholog_list;
  }
  my $helper_verb = "";
  if ($number_of_orthologs < 2) {
      $helper_verb = "\, is involved in ";
  }
  if ($number_of_orthologs >= 2) {
      $helper_verb = "\, are involved in ";
  }
#
my %elegans_process_elements_hash=();
foreach my $element_id (@elegans_array){
 if ($gene_class_hash{$element_id}){
  my @gene_process_array=();
  my $gene_process = $elegans_process_hash{$element_id};
  if ($gene_process){
   @gene_process_array=split(/\,/, $gene_process);
  }

   my @uniq_gene_process_array=uniq(sort(@gene_process_array));
   $elegans_process_elements_hash{$element_id} = \@uniq_gene_process_array;
 }
}
#
my $process_intersection_array_ref = array_common_elements(\%elegans_process_elements_hash);
my @process_intersection_array = @{$process_intersection_array_ref};
my @uniq_process_intersection_array = uniq(sort(@process_intersection_array));
my $gene_process = join(',', @uniq_process_intersection_array);
#
 if (@uniq_process_intersection_array > 0){
 my $gene_row = $gene_id . "\," . $gene_name . "\," . $gene_process;
# check for granularity
 my $processes_ref = ConciseDescriptions::get_granular_array($ontology, $gene_row, "IEA", "ISS", \%parents);
 my @processes = @{$processes_ref};
 my @uniq_processes = uniq(sort(@processes));
 my $granular_process = join(',', @uniq_processes);

 if (@uniq_processes){
# convert the ids into terms and put them in the array
    my $process_count = 0;
    foreach my $process (@uniq_processes){
     next if ($process =~/$gene_id/);
     next if ($process =~/$gene_name/);
     next if ($process =~/\[IEA\]/);
     next if ($process =~/\[ISS\]/);
     $process_count++;
     if ($process_count > 0){
        $process =~ s/^\s+//;
        $process =~ s/\s+$//;
     my $go_term = "";
             my ($evidence) = $process =~ /\[(.*)\]/;
             if ($evidence) {
                 $process =~ s/\[$evidence\]//g;
             }
             $process =~ s/^\s+//;
             $process =~ s/\s+$//;
           
     $go_term = ConciseDescriptions::get_ontology_term($process, \%gene_ontology);
     if ($go_term){
     if (lc($go_term) !~/obsolete/){
         $go_term .= "\[$evidence\]";
       if (grep {$_ =~/\Q$go_term\E/} @GO_processes) {
           print "$go_term ignored\; already in array\n";
        } else {
           push(@GO_processes, $go_term);
       }
      }
     }
    }
   }
     
     if (@GO_processes){
      my $go_process_print = join(',', @GO_processes);

#
# yet untested, $imp="IMP" would identify processes with evidence code separately; else it includes them without 
# adding information that it was "based on mutant phenotype"
#
     my $imp = "";
     my $process_goterm = get_verb_process_goterm(\@GO_processes, $imp);
     if ($process_goterm) {   
        $process_string = $process_goterm;
        $process_string =~ s/^\s+//;
        $process_string =~ s/\s+$//;
        $process_string =~ s/ \,/\,/g;
        $process_string =~ s/ +/ /g;
#
# Remove evidence codes
#
        my ($evidence_to_remove) = $process_string =~ /\[(.*)\]/;
        if ($evidence_to_remove) {
             $process_string =~ s/\[[^\]]*\]//g;
           }
         $process_string =~ s/$molting_0/$molting/g;
         $process_string =~ s/$molting_1/$molting/g;
         $process_string =~ s/$molting_2/$molting/g;
         $process_string =~ s/$molting_3/$molting/g;

         $process_string =~ s/$embryo_1/$embryo/g;

         $process_string =~ s/$the_the/$the/g;
         $process_string =~ s/$blank_comma/$comma/g;

        $process_string =~ s/$synaptic_1_0/$synaptic_1_1/g;
        $process_string =~ s/$synaptic_2_0/$synaptic_2_1/g;
        $process_string =~ s/$synaptic_3_0/$synaptic_3_1/g;
        $process_string =~ s/$synaptic_4_0/$synaptic_4_1/g;
        $process_string =~ s/$synaptic_5_0/$synaptic_5_1/g;
        $process_string =~ s/ +/ /g;
     } # process_goterm
    } #GO_processes
 } # gene_process
  my $sentence = "";
  if (@GO_processes){
    $sentence = $helper_verb . $process_string;
  } else {
    $sentence = "";
  }
          my $doublespace = "  ";
          my $space = " ";
          $sentence =~ s/$structural/$is_structural/g;
# remove redundant is
          $sentence =~s/$is_is/$is/g;
          $sentence =~s/$is_a_is_a/$is_a/g;
# remove any double spaces
          $sentence =~s/$doublespace/$space/g;
          $sentence =~s/$growth_0/$growth_1/g;
# add semi-colon at end
#          $sentence .= "\;";
      if (length($sentence) > 0){
          $sentence .= "\;";
          $sentence =~s/$the_cell_0/$the_cell_1/g;
       my $out = $individual_path . $gene_id;
       my $gene_sentence = $ortholog_sentence . " " . $sentence;
          $gene_sentence =~s/$doublespace/$space/g;
          $gene_sentence =~s/$space_comma/$just_comma/g;
        write_file($out, $gene_sentence);
        write_file($output_file, {append => 1 }, $gene_id);
        write_file($output_file, {append => 1 }, "\n");
        write_file($output_file, {append => 1 }, $gene_sentence);
        write_file($output_file, {append => 1 }, "\n\n\n");
       }
  }
 }
}
exit(0);
sub unique {
 my @unique=();
 my $array_ref = shift;
 if ($array_ref){
 my @array = @{$array_ref};
 my %hash   = map { $_ => 1 } @array;
    @unique = keys %hash;
 }
 return \@unique;
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
      return \@gene_array;
}
sub get_wbgene_id_elegans_hash{
 my $file = shift;
 my $db_gene_list = shift;
 my %hash=();
 my @lines = read_file($file);
  foreach my $line (@lines){
  chomp($line);
  my @fields = split(/\t/, $line);
  my $gene_name = $fields[2];
  my $gene_id = $fields[1];
  if ($gene_id){
       $hash{$gene_id} = $gene_name;
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
        $hash{$WBGene} = $Gene;
     }
 }
 return \%hash;
}
sub get_wbgene_name_elegans_hash{
 my $file = shift;
 my $db_gene_list = shift; 
 my %hash=();
 my @lines = read_file($file);
  foreach my $line (@lines){
  chomp($line);
  my @fields = split(/\t/, $line);
  my $gene_name = $fields[2];
  my $gene_id = $fields[1];
  if ($gene_name){
       $hash{$gene_name} = $gene_id;
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
  if ($Gene){
        $hash{$Gene} = $WBGene;
    }
     }
 }
 return \%hash;
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
      if (($WBGene) and ($Gene)) {    
        $WBGene =~ s/^\s+//;
        $WBGene =~ s/\s+$//;
        $Gene =~ s/^\s+//;
        $Gene =~ s/\s+$//;
        $gene_name_hash{$WBGene} = $Gene;
       }
     }
 }
 my @db_lines = read_file($db_gene_list);
 foreach my $db_line (@db_lines){
  chomp($db_line);
  if ($db_line =~/WBGene/){
     my ($WBGene, $Gene) = split(/\t/, $db_line);
      if (($WBGene) and ($Gene)) {   
        $WBGene =~ s/^\s+//;
        $WBGene =~ s/\s+$//;
        $Gene =~ s/^\s+//;
        $Gene =~ s/\s+$//;
        $gene_name_hash{$WBGene} = $Gene;
      }
     }
 }
      return \%gene_name_hash;
}
sub get_gene_go_hash{
 my $file = shift;
 my $go = shift;
 my $ec = shift;
 my %hash=();
 my @exp_ec = ("EXP", "IDA", "IPI", "IMP", "IGI", "IEP");

 $go =~ s/^\s+//;
 $go =~ s/\s+$//;
 my @lines = read_file($file);
  foreach my $line (@lines){
   next if ($line =~/\!/);
   next if ($line =~/UniProt/);
   next if ($line =~/IBA/);
   next if ($line =~/IBD/);
   next if ($line =~/GO\:0003674/); # veto molecular_function
   next if ($line =~/GO\:0005554/); # veto molecular_function alt_id
   next if ($line =~/GO\:0005575/); # veto cellular_component
   next if ($line =~/GO\:0008372/); # veto cellular_component alt_id
   next if ($line =~/GO\:0008150/); # veto biological_process
   next if ($line =~/GO\:0000004/); # veto biological_process alt_id
   next if ($line =~/GO\:0007582/); # veto biological_process alt_id
   next if ($line =~/GO\:0005488/); # veto binding
   next if ($line =~/GO\:0005515/); # veto protein binding

   chomp($line);
   my @fields = split(/\t/, $line);
   my $gene_id = $fields[1]; 
   my $go_term = $fields[4]; 
   my $term;
   $go_term =~ s/^\s+//;
   $go_term =~ s/\s+$//;
   $gene_id =~ s/^\s+//;
   $gene_id =~ s/\s+$//;
   my $evidence_code = $fields[6]; 
   my $with = $fields[7];
   my $GO = $fields[8];
   if (($with =~/WBPhenotype/) and ($evidence_code=~/IEA/) and ($GO =~/P/)) {
         $evidence_code = "IMP";
   }
   if ($ec){
       $term = $go_term . "\[" . $evidence_code . "\]";
   } else {
       $term = $go_term;
   }

   if ($go_term){
   if ($go =~/$GO/){
    if ($hash{$gene_id}){
     if ($hash{$gene_id}!~/$go_term/g){
       if (($evidence_code) and ($ec)) {
         $term = $go_term . "\[" . $evidence_code . "\]";
       } else {
         $term = $go_term;
       }
        $hash{$gene_id} .= "\, " . $term;
     } elsif ((grep {$_ =~ $evidence_code} @exp_ec) and ($ec)) {
         my $temp_string = $hash{$gene_id};
         my $replace = $go_term . "\[" . $evidence_code . "\]";
            $temp_string =~ s/$go_term\[[^\]]*\]/$replace/g;
            $hash{$gene_id} = $temp_string;
     }
    } else {
       if (($evidence_code) and ($ec)) {
         $term = $go_term . "\[" . $evidence_code . "\]";
       } else{
         $term = $go_term;
       }
        $hash{$gene_id} = $term;
    }
   }
  }

 } 
 return \%hash;
}
sub get_verb_process_goterm{
    my $goterms_ref = shift;
    my $imp = shift;
    my $imp_string = " based on mutant phenotype";
    my @goterms = uniq(@{$goterms_ref});
    my $verb_goterm="";
    my $go_terms_print = join(',', @goterms);
# 
# The following was added (untested) to incorporate an IMP, "based on mutant phenotype".
# cluster the GO terms

   if ($imp !~/IMP/){
        foreach my $element (@goterms){
          my ($evidence_to_remove) = $element =~ /\[(.*)\]/;
          if ($evidence_to_remove) {
             $element =~ s/\[[^\]]*\]//g;
           }
        }
       $verb_goterm = get_verb_goterm_process(\@goterms);
 } else {     
    my @IMP=();
    my @other_ec =();
    if (@goterms){
       foreach my $goterm (@goterms){
         if ($goterm =~/\[IMP\]/){
            my ($evidence_to_remove) = $goterm =~ /\[(.*)\]/;
            if ($evidence_to_remove) {
                $goterm =~ s/\[[^\]]*\]//g;
            }
            if (grep {$_ =~/\Q$goterm\E/} @IMP) {
             print "$goterm ignored\; already in array\n";
            } else {
             push(@IMP, $goterm);
            }

         } else {
            my ($evidence_to_remove) = $goterm =~ /\[(.*)\]/;
            if ($evidence_to_remove) {
                $goterm =~ s/\[[^\]]*\]//g;
             }
            if (grep {$_ =~/\Q$goterm\E/} @other_ec) {
                print "$goterm ignored\; already in array\n";
            } else {
                push(@other_ec, $goterm);
            }
         }
      } 
    }
 my $phrase_imp = "";
 if (@IMP){
   $phrase_imp = get_verb_goterm_process(\@IMP);
   $phrase_imp .= $imp_string;
 }
 my $phrase_other_ec="";
 if (@other_ec){
   $phrase_other_ec = get_verb_goterm_process(\@other_ec);
 }
 if (($phrase_imp) and ($phrase_other_ec)){
    $verb_goterm = $phrase_other_ec . " and " . $phrase_imp;
 } elsif ($phrase_imp){
    $verb_goterm = $phrase_imp;
 } elsif ($phrase_other_ec){
    $verb_goterm = $phrase_other_ec;
  }
 }
 return $verb_goterm;
}
sub get_verb_goterm_process{
    my $goterms_ref = shift;
    my @goterms = uniq(@{$goterms_ref});
    my $go_terms_print = join(',', @goterms);

    my $helper_string="";
    my $verb_goterm="";
       if (@goterms){
         my $count = -1;
         my $number_of_terms = @goterms;
        foreach my $goterm (@goterms){
         chomp($goterm);
         $goterm =~ s/^\s+//;
         $goterm =~ s/\s+$//;
         $count++;
        if ($number_of_terms == 1){
         $verb_goterm = $goterm;
        } 
        if ($number_of_terms == 2){ 
           if ($count > 0) {
             $verb_goterm .= " and ";
             $verb_goterm .= $goterm;
           } else {
             $verb_goterm = $goterm;
           }
        } 
        if ($number_of_terms >= 3) {
         if ($count == 0) {
             $verb_goterm = $goterm;
         } 
         if (($count > 0) and ($count < ($number_of_terms-1))) {
             $verb_goterm .= "\, ";
             $verb_goterm .= $goterm;
         } 
         if ($count == ($number_of_terms-1)) {
             $verb_goterm .= " and ";
             $verb_goterm .= $goterm;
            }
         }
        } # goterms loop
       } # goterms gt 0
 return $verb_goterm; 
}
sub get_verb_function_goterm{
       my $goterms_ref = shift;
       my $ec = shift;

       my $helper_string_before_iea = " is predicted to have ";
       my $helper_string_before_iea_1 = " is predicted to be ";
       my $helper_string_after_iea = "\, based on protein domain information";

       my $helper_string_before_iss = " is predicted to have ";
       my $helper_string_after_iss = "\, based on sequence information";

       my $helper_string_before_exp = " exhibits ";
#      as per R. Kishore's changes
#       my $helper_string_after_exp = "\, based on experimental evidence";
       my $helper_string_after_exp = "";

       my $helper_string_before=""; 
       my $helper_string_after=""; 

       if ($ec =~ /IEA/){
        $helper_string_before=$helper_string_before_iea; 
        $helper_string_after=$helper_string_after_iea;
       } 
       if ($ec =~ /ISS/){
        $helper_string_before=$helper_string_before_iss; 
        $helper_string_after=$helper_string_after_iss;
       } 
       if ($ec =~ /EXP/){
        $helper_string_before=$helper_string_before_exp; 
        $helper_string_after=$helper_string_after_exp;
       } 

       my $verb_goterm="";
       my $watch_string_0 = "activity";
       my $watch_string_1 = "structural constituent";

       my $helper_string_0=$helper_string_before;
       my $helper_string_1=$helper_string_after;

       my @goterms = uniq(@{$goterms_ref});
       if (@goterms){
         my $count = 0;
         my $verb_count = 0;
         my $goterm="";
        foreach $goterm (@goterms){
        chomp($goterm);
        if (@goterms == 1){

         $helper_string_0=$helper_string_before;
         $helper_string_1=$helper_string_after;

           if ( $goterm =~/$watch_string_1/){
                $helper_string_0 = " is a ";
                $helper_string_1 = $helper_string_after;
           }

         $verb_count++;
         $verb_goterm =  $helper_string_0 . $goterm . $helper_string_1;      
       } elsif (@goterms == 2){ 

         $helper_string_0=$helper_string_before;
         $helper_string_1=$helper_string_after;

           if ( $goterm !~/$watch_string_0/){
                $helper_string_0 = "";
                $helper_string_1 = "";
            }

           if ($count > 0) {

           if ( $goterm !~/$watch_string_0/){
                $helper_string_0 = "";
                $helper_string_1 = "";
            }

            if ($verb_goterm =~/$helper_string_0/){
             $verb_goterm .= " and ";
             $verb_goterm .= $goterm;
             $verb_goterm .= $helper_string_1; 
            } else {
             $verb_goterm .= " and ";
             $verb_goterm .= $helper_string_0;    
             $verb_goterm .= $goterm;        
             $verb_goterm .= $helper_string_1; 
            }

           } else {

           if ( $goterm !~/$watch_string_0/){
                $helper_string_0 = "";
            }

           if ( $goterm =~/$watch_string_1/){
                $helper_string_0 = "";
                $helper_string_1 = $helper_string_after;
           }

               $verb_goterm = $helper_string_0 . $goterm ;
           }

       } elsif (@goterms > 2) {

           if ( $goterm !~/$watch_string_0/){
                $helper_string_0 = "";
                $helper_string_1 = "";
            }


         if ($count == 0) {
            if ( $goterm =~/$watch_string_1/){
                $helper_string_0 = "";
                $helper_string_1 = " and " . $helper_string_before;
           }

           if ( $goterm =~/$watch_string_1/){
                            $verb_goterm = $helper_string_0 . $goterm . $helper_string_1;
           } else {
                             $verb_goterm = $helper_string_0 . $goterm . "\, ";
           }
         } elsif ($count < (@goterms-1)){

             if (($verb_goterm =~/$helper_string_1/) and ($goterm !~/$watch_string_1/)) {
                $verb_goterm .= $goterm . "\, ";
            } elsif ($verb_goterm =~/$helper_string_0/) {
             $verb_goterm .= "\, ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= "\, and "; 
             $verb_goterm .= $helper_string_0;    
             $verb_goterm .= $goterm;             
             $verb_goterm .= $helper_string_1;  
            }
         } elsif ($count >= (@goterms-1)) {
                  $helper_string_0=$helper_string_before;
                  $helper_string_1=$helper_string_after;

            if (($verb_goterm =~/$helper_string_0/) and ($goterm !~/$watch_string_1/)) {
                 $verb_goterm .= "\, and ";
                 $verb_goterm .= $goterm;
                 $verb_goterm .= $helper_string_1; 
            } elsif ($verb_goterm =~/$helper_string_0/){
             $verb_goterm .= "\, and ";
             $verb_goterm .= $goterm;
             $verb_goterm .= $helper_string_1; 
            } elsif (($verb_goterm !~/$helper_string_1/) or ($goterm =~/$watch_string_1/)) {
             $verb_goterm .= "\, and ";
             $verb_goterm .= $helper_string_0;    
             $verb_goterm .= $goterm;             
             $verb_goterm .= $helper_string_1;  
            }

         } 
        }
         $count++;
        } # goterms loop
       } # goterms gt 0 
 return $verb_goterm;
}
sub get_verb_component_goterm{
       my $goterms_ref = shift;
       my @goterms = uniq(@{$goterms_ref});
       my $verb_goterm="";
       my $helper_string = "";
       my $watch_string = "integral component of";
       my $watch_string_the = "integral component of the";
       my $helper_string_0 = " is localized to the ";
       my $helper_string_1 = " is an ";
       if (@goterms){
         my $count = 0;
         my $verb_count = 0;
         my $goterm="";
        foreach $goterm (@goterms){
        if (@goterms == 1){
         $helper_string = $helper_string_0;
           if ( $goterm =~/$watch_string/){
                $helper_string = $helper_string_1;
           }
         $verb_count++;
         $verb_goterm =  $helper_string . $goterm;      
       } elsif (@goterms == 2){ 

           $helper_string = $helper_string_0;
           if ( $goterm =~/$watch_string/){
                $helper_string = $helper_string_1;
           }
           if ($count > 0) {
            if ($verb_goterm =~/$helper_string/){
             $verb_goterm .= " and the ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= " and ";
             $verb_goterm .= $helper_string;    
             $verb_goterm .= $goterm;               
            }

           } else {
             $verb_goterm = $helper_string . $goterm;
           }

       } elsif (@goterms > 2) {
            $helper_string = $helper_string_0;
           if ( $goterm =~/$watch_string/){
                $helper_string = $helper_string_1;
           }
         if ($count == 0) {
             $verb_goterm = $helper_string . $goterm;
         } elsif ($count < (@goterms-1)) {

            if ($verb_goterm =~/$helper_string/){
             $verb_goterm .= "\, the ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= " and "; 
             $verb_goterm .= $helper_string;    
             $verb_goterm .= $goterm;               
            }

         } elsif ($count >= (@goterms-1)) {
            if ($verb_goterm =~/$helper_string/){
             $verb_goterm .= " and the ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= " and ";
             $verb_goterm .= $helper_string;    
             $verb_goterm .= $goterm;               
            }
         } 
        }
         $count++;
        } # goterms loop
       } # goterms gt 0 
 return $verb_goterm;
}
#
sub get_exp_go_array{
 my $go_array_ref = shift;
 my @array=();
 my @go_array = uniq(@{$go_array_ref});
  foreach my $go (@go_array){
        $go =~ s/^\s+//;
        $go =~ s/\s+$//;
   if ( ($go =~/\[EXP\]/) or 
  ($go =~/\[IDA\]/) or ($go =~/\[IPI\]/) or 
  ($go =~/\[IMP\]/) or ($go =~/\[IMP\]/) or 
  ($go =~/\[IGI\]/) or ($go =~/\[IEP\]/) ) {
    if (grep {$_ =~/\Q$go\E/} @array) {
        print "$go ignored\; already in array\n";
      } else {
        push(@array, $go);
      }  
  }
 }

  my @sorted_array = uniq(sort(@array));
  
 return \@sorted_array;
}
sub get_stat_go_array{
 my $go_array_ref = shift;
 my @array=();
 my @go_array = uniq(@{$go_array_ref});
  foreach my $go (@go_array){
        $go =~ s/^\s+//;
        $go =~ s/\s+$//;
   if ( ($go =~/\[ISS\]/) or 
  ($go =~/\[ISO\]/) or ($go =~/\[ISA\]/) or 
  ($go =~/\[ISM\]/) or ($go =~/\[IGC\]/) or 
  ($go =~/\[IBA\]/) or ($go =~/\[IBD\]/) or 
  ($go =~/\[IKR\]/) or ($go =~/\[IRD\]/) or
  ($go =~/\[RCA\]/) ) {
    if (grep {$_ =~/\Q$go\E/} @array) {
        print "$go ignored\; already in array\n";
      } else {
        push(@array, $go);
      }  
  }
 }
  my @sorted_array = uniq(sort(@array));
  
 return \@sorted_array;
}
sub get_iea_go_array{
 my $go_array_ref = shift;
 my @array=();
 my @go_array = uniq(@{$go_array_ref});
  foreach my $go (@go_array){
        $go =~ s/^\s+//;
        $go =~ s/\s+$//;
   if (($go =~/\[IEA\]/) or 
  ($go =~/\[TAS\]/) or ($go =~/\[NAS\]/) or 
  ($go =~/\[IC\]/) or ($go =~/\[ND\]/)) {
    if (grep {$_ =~/\Q$go\E/} @array) {
        print "$go ignored\; already in array\n";
      } else {
        push(@array, $go);
      }
  }
 }
#  my $sorted_array_ref = unique(\@array);
#  my @sorted_array = sort(@{$sorted_array_ref});
   my @sorted_array = uniq( sort(@array) );
 
 return \@sorted_array;
}
sub get_function_string {
 my $functions_array_ref = shift;
 my $ec = shift;
 my $gene_ontology_ref = shift;
 my $function_string = "";
 my @functions = uniq(@{$functions_array_ref});
 my @GO_functions = ();
 
 if (@functions){
    my $function_count=0;
    foreach my $function (@functions){
          next if ($function =~/GO\:0005488/); # ignore binding
          next if ($function =~/GO\:0005515/); # ignore protein binding
          $function_count++;
          if ($function_count > 0){
          my ($evidence) = $function =~ /\[(.*)\]/;
             if ($evidence) {
                 $function =~ s/\[$evidence\]//g;
             }
             $function =~ s/^\s+//;
             $function =~ s/\s+$//;
     my $go_term = ConciseDescriptions::get_ontology_term($function, $gene_ontology_ref);  
     $go_term .= "[$evidence]";
     if ( $go_term !~ /activity/) { 
          $go_term =~s/binding/binding activity/g;
     }  
     if (grep {$_ =~/\Q$go_term\E/} @GO_functions) {
        print "$go_term ignored\; already in array\n";
      } else {
        push(@GO_functions, $go_term);
      }
    }
    } # foreach functions

    if (@GO_functions){
       my $gof_string = "";
       foreach my $gof (@GO_functions){
        $gof_string .= "$gof\t";
       }
       $function_string = get_verb_function_goterm(\@GO_functions, $ec);
      }
my @evidences=("EXP","IMP","IGI","IPI","IDA","IEP","IEA","ISS","ISA","ISO","ISM","IGC","IBA","IBD","IKR","IRD","RCA","TAS","NAS","ND","IC");
# Remove evidence codes
         foreach my $ec (@evidences){
	    $function_string =~ s/\[$ec\]//g;
         }
      $function_string =~ s/^\s+//;
      $function_string =~ s/\s+$//;
      $function_string =~ s/ \,/\,/g;
# Remove evidence codes
 my $size = @functions;

 } # if gene_functions array

 return $function_string;
}
sub array_common_elements{
my $array_ref_hash_ref = shift;
#
# use this routine like the following:
# my $intersection_ref = array_common_elements(\%array_ref_hash);
# my @intersection = @{$intersection_ref};
#
my @intersection = ();
my %records = %{$array_ref_hash_ref};
# records should be hash; 
# each key mapped value should be a reference to an array
my %count;
foreach my $arr_ref (values %records) {
    foreach my $elem (@$arr_ref) {
        $elem =~ s/^\s+//;
        $elem =~ s/\s+$//;
        $count{$elem}++;
    }
}

my $num_arrays = scalar(keys %records);
foreach my $elem (keys %count) {
        $elem =~ s/^\s+//;
        $elem =~ s/\s+$//;
    #If all the arrays contained this element, 
    #allowing for multiple entries per array
    if ($count{$elem} >= $num_arrays) {
      if (grep {$_ =~/\Q$elem\E/} @intersection) {
        print "$elem ignored\; already in array\n";
      } else {
        push(@intersection, $elem);
      }
    }
}
# returns a reference to an array of common values
 my @uniq_intersection = uniq(sort(@intersection));
return \@uniq_intersection;
}
sub get_id_name_hash{
 my $orthology=shift;
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
               chomp($a_line);
               $a_line =~ s/^\s+//;
               $a_line =~ s/\s+$//;

               my ($a_wb_gene_id, $a_gene_species_name, $a_id, $a_name, $a_source) = split(/\t/, $a_line);
               $a_id =~ s/^\s+//;
               $a_id =~ s/\s+$//;
               $a_name =~ s/^\s+//;
               $a_name =~ s/\s+$//;
               $a_wb_gene_id =~ s/^\s+//;
               $a_wb_gene_id =~ s/\s+$//;

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
               chomp($a_line);
               $a_line =~ s/^\s+//;
               $a_line =~ s/\s+$//;

               my ($a_wb_gene_id, $a_gene_species_name, $a_id, $a_name, $a_source) = split(/\t/, $a_line);
               my ($a_gene_name, $a_species_name) = split(/ /, $a_gene_species_name);
               $a_id =~ s/^\s+//;
               $a_id =~ s/\s+$//;
               $a_name =~ s/^\s+//;
               $a_name =~ s/\s+$//;
               $a_wb_gene_id =~ s/^\s+//;
               $a_wb_gene_id =~ s/\s+$//;

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
 return \%gene_id_hash, \%gene_name_hash;
}
