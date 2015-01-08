#!/usr/bin/env perl
use ConciseDescriptions;
use LWP::Simple;
use List::MoreUtils qw(uniq);
use List::Util qw/first/;
use File::Slurp;
use Text::CSV;
use warnings;
use strict;

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
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
my $gene_list_dir = $home . "release/$PRODUCTION_RELEASE/c_elegans/gene_lists/";
# Define GO term file
my $ontology = "GO";
my $go_file = $gene_list_dir . "newterms.txt";
my $go_altid =  $gene_list_dir . "goid_altid.txt";
my $input_ace_file = $gene_list_dir . "go_terms.ace";
#
# Create a file containing the GO terms and their respective parents
# build parent/child hash for granularity
my ($parents_ref, $children_ref) =  ConciseDescriptions::get_ontology_parents_children($input_ace_file, $ontology);
my %parents = %$parents_ref;
my %children = %$children_ref;
# Create GO hash so that terms can be referenced by GO ID.
my $gene_ontology_ref = ConciseDescriptions::get_ontology_hash($ontology, $go_file, $go_altid);
my %gene_ontology = %$gene_ontology_ref;
my $go_dir = $home . "release/$PRODUCTION_RELEASE/$species/gene_ontology/input_files/";
my $go_elegans_dir = $home . "release/$PRODUCTION_RELEASE/c_elegans/gene_ontology/input_files/";
my $gene_association_species_file = $go_dir . "gene_association.$RELEASE.wb.$species";
my $gene_association_elegans_file = $go_elegans_dir . "gene_association.$RELEASE.wb.c_elegans";
my $orthology = $home . "release/$PRODUCTION_RELEASE/$species/orthology/";
my $elegans_orthology = $home . "release/$PRODUCTION_RELEASE/c_elegans/orthology/";

my $go_path = $home . "release/$PRODUCTION_RELEASE/$species/gene_ontology/output_files/";
my $output_file = $go_path . "sentences_for_GO.txt";
my $xrefs = $orthology . "input_files/$species.$RELEASE.xrefs.txt";
my $xrefs_elegans = $elegans_orthology . "input_files/c_elegans.$RELEASE.xrefs.txt";
my $blastp = $orthology . "input_files/$species.$RELEASE.best_blastp_hits.txt";
#
my $species_ce_protein_hash_ref = get_species_ce_protein_blastp_hash($species, $blastp);
my %species_ce_protein_hash =%$species_ce_protein_hash_ref;
my $elegans_gene_name_xrefs_hash_ref = get_gene_name_xrefs_hash("c_elegans", $xrefs_elegans);
my %elegans_gene_name_xrefs_hash = %$elegans_gene_name_xrefs_hash_ref;
my $elegans_gene_protein_xrefs_hash_ref = get_gene_protein_xrefs_hash("c_elegans", $xrefs_elegans);
my %elegans_gene_protein_xrefs_hash = %$elegans_gene_protein_xrefs_hash_ref;
my $elegans_name_xrefs_hash_ref = get_inverse_gene_protein_xrefs_hash("c_elegans", $xrefs_elegans);
my %elegans_name_xrefs_hash = %$elegans_name_xrefs_hash_ref;
#
my $gene_protein_hash_ref = get_gene_protein_xrefs_hash($species, $xrefs); 
my %gene_protein_hash = %$gene_protein_hash_ref;
my $wbgene_elegans_name_hash_ref = get_wbgene_name_elegans_hash($gene_association_elegans_file);
my %wbgene_elegans_name_hash = %$wbgene_elegans_name_hash_ref;
my $wbgene_id_elegans_hash_ref = get_wbgeneid_elegans_hash($gene_association_elegans_file);
my %wbgene_id_elegans_hash = %$wbgene_elegans_name_hash_ref;
#
my $doublecomma = "\,\,";
my $the_the = "the the";
my $the = "the";
my $blank_comma = "\, \,";
my $comma = "\, ";
my $nadp =" NAD or NADP as acceptor\,";
my $component_watch_string = "integral component of";
my $component_watch_string_the = "integral component of the";
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
#
my $helper_string_after_iea = "\, based on protein domain information";
my $helper_string_before_iea = "based on protein domain information\, ";
#
if (-e $go_path){
   print "$go_path exists\n"; 
} else {
   mkdir $go_path or die "could not create $go_path";
}
my $individual_path = $go_path . "individual_gene_sentences/";
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

my $gene_array_ref = get_gene_array($gene_association_species_file);
my @gene_array = @$gene_array_ref;
my $gene_name_hash_ref = get_gene_name_hash($gene_association_species_file);
my %gene_name_hash = %$gene_name_hash_ref;
my $go = "P";
my $gene_process_hash_ref = get_gene_go_hash($gene_association_species_file, $go);
my %gene_process_hash = %$gene_process_hash_ref;
my $elegans_process_hash_ref = get_gene_go_hash($gene_association_elegans_file, $go, 1);
my %elegans_process_hash = %$elegans_process_hash_ref;
   $go = "C";
my $gene_component_hash_ref = get_gene_go_hash($gene_association_species_file, $go);
my %gene_component_hash = %$gene_component_hash_ref;
my $elegans_component_hash_ref = get_gene_go_hash($gene_association_elegans_file, $go, 1);
my %elegans_component_hash = %$elegans_component_hash_ref;
   $go = "F";
my $gene_function_hash_ref = get_gene_go_hash($gene_association_species_file, $go);
my %gene_function_hash = %$gene_function_hash_ref;
my $elegans_function_hash_ref = get_gene_go_hash($gene_association_elegans_file, $go, 1);
my %elegans_function_hash = %$elegans_function_hash_ref;

foreach my $gene_id (@gene_array){
 my @gene_processes = ();
 my @gene_functions = ();
 my @gene_components = ();

 my @GO_processes  = ();
 my @GO_functions  = ();
 my @GO_components = ();
 
 my $function_string = "";
 my $component_string = "";
 my $process_string = "";

 my $gene_name = $gene_name_hash{$gene_id};
 my $gene_process = $gene_process_hash{$gene_id};
 if ($gene_process){
 my $gene_row = $gene_id . "\," . $gene_name . "\," . $gene_process;
# print "gene row is $gene_row\n";
# check for granularity
 my $processes_ref = ConciseDescriptions::get_granular_array($ontology, $gene_row, "", "", \%parents);
 my @processes = @{$processes_ref};
 if (@processes){
# convert the ids into terms and put them in the array
    my $process_count = 0;
    foreach my $process (@processes){
     $process_count++;
     if ($process_count gt 2){
        $process =~ s/^\s+//;
        $process =~ s/\s+$//;
     my $go_term = ConciseDescriptions::get_ontology_term($process, \%gene_ontology);    
     push(@GO_processes, $go_term);
#     print "process go_term $process\t$go_term\n";
     }
    }
# 
#     my $process_goterm = join(', ', @GO_processes)
#
     if (@GO_processes){
     my $process_goterm = get_verb_process_goterm(\@GO_processes);
     if ($process_goterm) {   
        $process_string = $process_goterm;
#        print "process string = $process_string\n";
#        $process_string .= "\;\n";
        $process_string =~ s/^\s+//;
        $process_string =~ s/\s+$//;
        $process_string =~ s/ \,/\,/g;
        $process_string =~ s/ +/ /g;
# Remove evidence codes
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
     }# process_goterm
    }
 } # gene_process
}
 my $gene_function = $gene_function_hash{$gene_id};
 if ($gene_function){
 my $gene_row = $gene_id . "\," . $gene_name . "\," . $gene_function;
# check for granularity
 my $functions_ref = ConciseDescriptions::get_granular_array($ontology, $gene_row, "", "", \%parents);
 my @functions = @{$functions_ref};
 if (@functions){
    my $function_count=0;
    foreach my $function (@functions){
          $function_count++;
          if ($function_count gt 2){
          my ($evidence) = $function =~ /\[(.*)\]/;
             if ($evidence) {
                 $function =~ s/\[$evidence\]//g;
             }
             $function =~ s/^\s+//;
             $function =~ s/\s+$//;
     my $go_term = ConciseDescriptions::get_ontology_term($function, \%gene_ontology);  
     if ( $go_term !~ /activity/) { 
          $go_term =~s/binding/binding activity/g;
     }  
     push(@GO_functions, $go_term);
    }
    } # foreach functions
    if (@GO_functions){
      $function_string = get_verb_function_goterm(\@GO_functions);
#      $function_string .= "\;\n";
      $function_string =~ s/^\s+//;
      $function_string =~ s/\s+$//;
      $function_string =~ s/ \,/\,/g;
      $function_string =~ s/ +/ /g;
      $function_string =~ s/$the_the/$the/g;
      $function_string =~ s/$blank_comma/$comma/g;
      $function_string =~ s/$nadp//g;
      $function_string =~ s/$doublecomma/$comma/g;     
      $function_string =~ s/ +/ /g;
     }
# Remove evidence codes
 } # if gene_functions array
}
 my $gene_component = $gene_component_hash{$gene_id};
 if ($gene_component){
  my $gene_row = $gene_id . "\," . $gene_name . "\," . $gene_component;
 if ($gene_component){
    @gene_components = split(/\t/, $gene_component);
 }
   my $components_ref = ConciseDescriptions::get_granular_array($ontology, $gene_row, "", "", \%parents);
   my @components = @{$components_ref};
   my $component_count = 0;
   if (@components) {
    foreach my $component (@components){
    $component_count++;
    if ($component_count gt 2){
             my ($evidence) = $component =~ /\[(.*)\]/;
             if ($evidence) {
                 $component =~ s/\[$evidence\]//g;
             }
             $component =~ s/^\s+//;
             $component =~ s/\s+$//;
           
     my $go_term = ConciseDescriptions::get_ontology_term($component, \%gene_ontology);  
     push(@GO_components, $go_term);
    }
   } # foreach components
       if (@GO_components){
        $component_string = get_verb_component_goterm(\@GO_components);

        $component_string =~ s/$component_watch_string/$component_watch_string_the/g;
        $component_string =~ s/$the_the/$the/g;
        $component_string =~ s/$blank_comma/$comma/g;

        $component_string =~ s/ +/ /g;
#        $component_string .= "\n";
        $component_string =~ s/^\s+//;
        $component_string =~ s/\s+$//;
        $component_string =~ s/ \,/\,/g;
        $component_string =~ s/ +/ /g;
      }
   } # if components array
 }
  my $sentence = "";
  if ((@GO_components) and (@GO_processes) and (@GO_functions)){
    $sentence = $gene_name . " " . $process_string . "\, " . $function_string . " and " . $component_string;
  } elsif ((@GO_processes) and (@GO_functions)){
    $sentence = $gene_name . " " . $process_string . " and " . $function_string;
  } elsif ((@GO_components) and (@GO_functions)){
    $sentence = $gene_name . " " . $function_string . " and " . $component_string;
  } elsif ((@GO_components) and (@GO_processes)){
    $sentence = $gene_name . " " . $process_string . " and " . $component_string;
  } elsif (@GO_components){
    $sentence = $gene_name . " " . $component_string;
  } elsif (@GO_processes){
    $sentence = $gene_name . " " . $process_string;
  } elsif (@GO_functions){
    $sentence = $gene_name . " " . $function_string;
  } else {
    $sentence = "";
  }
# add semi-colon at the end
#    $sentence .= "\;";
#  print "$sentence\n\n";
          my $doublespace = "  ";
          my $space = " ";
# remove any double spaces
          $sentence =~s/$doublespace/$space/g;
# if sentence has data based on protein domain information, that information should be first.
         if ($species !~/elegans/){
          if ($sentence =~/$helper_string_after_iea/){
#           print "old sentence\: $sentence\n";
           my $tmp = $sentence;
           $tmp =~s/$helper_string_after_iea//g;
           $sentence = $helper_string_before_iea . $tmp;
#           print "new sentence\: $sentence\n";
          }
         }
# add semi-colon at end
#          $sentence .= "\;";
# add ortholog sentence
          my $ortholog_sentence = get_ortholog_sentence($gene_id, $species_prefix);
#          print "sentence is $gene_id\t$gene_name\t$ortholog_sentence\n";
      if (length($sentence) gt 0){
          if (($ortholog_sentence ne "") and ($species !~/elegans/)){
              $sentence .= "\; in C\. elegans\, " . $ortholog_sentence;
           } 
              $sentence .= "\.";
      
       my $out = $individual_path . $gene_id;
        write_file($out, $sentence);
        write_file($output_file, {append => 1 }, $sentence);
        write_file($output_file, {append => 1 }, "\n\n\n");
       }
 }
exit(0);

sub get_wbgeneid_elegans_hash{
 my $file = shift;
 my %hash=();
 my @lines = read_file($file);
  foreach my $line (@lines){
  my @fields = split(/\t/, $line);
  my $elegans_name = $fields[4];
  my $gene_id = $fields[0];
  if ($elegans_name){
   if ($hash{$gene_id}){
       $hash{$gene_id} .= "\," . $elegans_name;
   } else {
       $hash{$gene_id} = $elegans_name;
   }
  }
 }
 return \%hash;
}
sub get_wbgene_name_elegans_hash{
 my $file = shift;
 my %hash=();
 my @lines = read_file($file);
  foreach my $line (@lines){
  my @fields = split(/\t/, $line);
  my $gene_name = $fields[2];
  my $gene_id = $fields[1];
  if ($gene_name){
       $hash{$gene_name} = $gene_id;
   }
 }
 return \%hash;
}
sub get_gene_array{
 my $file = shift;
 my @array=();
 my @lines = read_file($file);
  foreach my $line (@lines){
   next if ($line =~/\!/);
   my @fields = split(/\t/, $line);
   my $gene_id = $fields[1];
      $gene_id =~ s/^\s+//;
      $gene_id =~ s/\s+$//;
   if ($gene_id){
    push(@array, $gene_id);  
  }
 }
  my @sorted_array = sort(uniq(@array));
  
 return \@sorted_array;
}
sub get_gene_name_hash{
 my $file = shift;
 my %hash=();
 my @lines = read_file($file);
  foreach my $line (@lines){
   next if ($line =~/\!/);
   my @fields = split(/\t/, $line);
   my $gene_id = $fields[1];
   my $gene_name = $fields[2];
      $gene_name =~ s/^\s+//;
      $gene_name =~ s/\s+$//;
      $gene_id   =~ s/^\s+//;
      $gene_id   =~ s/\s+$//;
   if ($gene_id){
    $hash{$gene_id} = $gene_name;
  }
 }
 return \%hash;
}
sub get_gene_go_hash{
 my $file = shift;
 my $go = shift;
 my $ec = shift;
 my %hash=();
 
 $go =~ s/^\s+//;
 $go =~ s/\s+$//;
 my @lines = read_file($file);
  foreach my $line (@lines){
   next if ($line =~/\!/);
   my @fields = split(/\t/, $line);
   my $gene_id = $fields[1]; 
   my $go_term = $fields[4]; 
   my $term;
   $go_term =~ s/^\s+//;
   $go_term =~ s/\s+$//;
   $gene_id =~ s/^\s+//;
   $gene_id =~ s/\s+$//;
   my $evidence_code = $fields[6]; 
   my $GO = $fields[8]; 
   if ($ec){
       $term = $go_term . "\[" . $evidence_code . "\]";
   } else{
       $term = $go_term;
   }
   if ($term){
    if (($term !~/GO\:0005488/)and($term !~/GO\:0005515/)){
   if ($go =~/$GO/){
    if ($hash{$gene_id}){
        $hash{$gene_id} .= "\, " . $term;
    } else {
        $hash{$gene_id} = $term;
    }
   }
  }
  }
 } 
 return \%hash;
}

sub get_verb_process_goterm{
    my $goterms_ref = shift;
    my $helper_string_0 = " is involved in ";
    my $helper_string_1 = " functions in ";
    my $watch_string = "involved in";

    my @goterms = @$goterms_ref;
    my $helper_string="";
    my $verb_goterm="";
       if (@goterms){
         my $count = 0;
         my $verb_count = 0;
        foreach my $goterm (@goterms){
#        chomp($goterm);
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
           if ($count gt 0) {
            if ($verb_goterm =~ $helper_string){
             $verb_goterm .= " and ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= " and ";
             $verb_goterm .= $helper_string;    
             $verb_goterm .= $goterm;               
            }

           } else {
             $verb_goterm = $helper_string . $goterm;
           }

       } elsif (@goterms gt 2) {
            $helper_string = $helper_string_0;
           if ( $goterm =~/$watch_string/){
                $helper_string = $helper_string_1;
           }
         if ($count eq 0) {
             $verb_goterm = $helper_string . $goterm;
         } elsif ($count lt  @goterms-1) {

            if ($verb_goterm =~ $helper_string){
             $verb_goterm .= "\, ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= " and "; 
             $verb_goterm .= $helper_string;    
             $verb_goterm .= $goterm;               
            }

         } elsif ($count ge @goterms-1) {
            if ($verb_goterm =~ $helper_string){
             $verb_goterm .= " and ";
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
sub get_verb_function_goterm_exp{
       my $goterms_ref = shift;
#
       my $helper_string_before_exp = " exhibits ";
#      as per R. Kishore's changes
#       my $helper_string_after_exp = "\, based on experimental evidence";
       my $helper_string_after_exp = "";
       my $helper_string_before=$helper_string_before_exp; 
       my $helper_string_after=$helper_string_after_exp; 

#       my $helper_string_before = shift;
#       my $helper_string_after  = shift;
       my $verb_goterm="";
       my $watch_string_0 = "activity";
       my $watch_string_1 = "structural constituent";

       my $helper_string_0=$helper_string_before;
       my $helper_string_1=$helper_string_after;

       my @goterms = @{$goterms_ref};
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

           if ($count gt 0) {

           if ( $goterm !~/$watch_string_0/){
                $helper_string_0 = "";
                $helper_string_1 = "";
            }

            if ($verb_goterm =~ $helper_string_0){
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

       } elsif (@goterms gt 2) {

           if ( $goterm !~/$watch_string_0/){
                $helper_string_0 = "";
                $helper_string_1 = "";
            }


         if ($count eq 0) {
            if ( $goterm =~/$watch_string_1/){
                $helper_string_0 = "";
                $helper_string_1 = " and " . $helper_string_before;
           }

           if ( $goterm =~/$watch_string_1/){
                            $verb_goterm = $helper_string_0 . $goterm . $helper_string_1;
           } else {
                             $verb_goterm = $helper_string_0 . $goterm . "\, ";
           }
         } elsif ($count lt  @goterms-1){

             if (($verb_goterm =~ $helper_string_1) and ($goterm !~ $watch_string_1)) {
                $verb_goterm .= $goterm . "\, ";
            } elsif ($verb_goterm =~ $helper_string_0) {
             $verb_goterm .= "\, ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= "\, and "; 
             $verb_goterm .= $helper_string_0;    
             $verb_goterm .= $goterm;             
             $verb_goterm .= $helper_string_1;  
            }
         } elsif ($count ge @goterms-1) {
                  $helper_string_0=$helper_string_before;
                  $helper_string_1=$helper_string_after;

            if (($verb_goterm =~ $helper_string_0) and ($goterm !~ $watch_string_1)) {
                 $verb_goterm .= "\, and ";
                 $verb_goterm .= $goterm;
                 $verb_goterm .= $helper_string_1; 
            } elsif ($verb_goterm =~ $helper_string_0){
             $verb_goterm .= "\, and ";
             $verb_goterm .= $goterm;
             $verb_goterm .= $helper_string_1; 
            } elsif (($verb_goterm !~ $helper_string_1) or ($goterm =~ $watch_string_1)) {
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
sub get_verb_function_goterm{
       my $goterms_ref = shift;
# assuming that the evidence codes are IEA for now.
       my $helper_string_before_iea = " is predicted to have ";
       my $helper_string_before_iea_1 = " is predicted to be ";
       my $helper_string_before_iss = " is predicted to have ";
       my $helper_string_after_iea = "\, based on protein domain information";
       my $helper_string_after_iss = "\, based on sequence information";

       my $helper_string_before=$helper_string_before_iea; 
       my $helper_string_after=$helper_string_after_iea; 

#       my $helper_string_before = shift;
#       my $helper_string_after  = shift;
       my $verb_goterm="";
       my $watch_string_0 = "activity";
       my $watch_string_1 = "structural constituent";

       my $helper_string_0=$helper_string_before;
       my $helper_string_1=$helper_string_after;

       my @goterms = @{$goterms_ref};
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

           if ($count gt 0) {

           if ( $goterm !~/$watch_string_0/){
                $helper_string_0 = "";
                $helper_string_1 = "";
            }

            if ($verb_goterm =~ $helper_string_0){
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

       } elsif (@goterms gt 2) {

           if ( $goterm !~/$watch_string_0/){
                $helper_string_0 = "";
                $helper_string_1 = "";
            }


         if ($count eq 0) {
            if ( $goterm =~/$watch_string_1/){
                $helper_string_0 = "";
                $helper_string_1 = " and " . $helper_string_before;
           }

           if ( $goterm =~/$watch_string_1/){
                            $verb_goterm = $helper_string_0 . $goterm . $helper_string_1;
           } else {
                             $verb_goterm = $helper_string_0 . $goterm . "\, ";
           }
         } elsif ($count lt  @goterms-1){

             if (($verb_goterm =~ $helper_string_1) and ($goterm !~ $watch_string_1)) {
                $verb_goterm .= $goterm . "\, ";
            } elsif ($verb_goterm =~ $helper_string_0) {
             $verb_goterm .= "\, ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= "\, and "; 
             $verb_goterm .= $helper_string_0;    
             $verb_goterm .= $goterm;             
             $verb_goterm .= $helper_string_1;  
            }
         } elsif ($count ge @goterms-1) {
                  $helper_string_0=$helper_string_before;
                  $helper_string_1=$helper_string_after;

            if (($verb_goterm =~ $helper_string_0) and ($goterm !~ $watch_string_1)) {
                 $verb_goterm .= "\, and ";
                 $verb_goterm .= $goterm;
                 $verb_goterm .= $helper_string_1; 
            } elsif ($verb_goterm =~ $helper_string_0){
             $verb_goterm .= "\, and ";
             $verb_goterm .= $goterm;
             $verb_goterm .= $helper_string_1; 
            } elsif (($verb_goterm !~ $helper_string_1) or ($goterm =~ $watch_string_1)) {
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
       my @goterms = @$goterms_ref;
       my $verb_goterm="";
       my $helper_string = "";
# assuming that the evidence codes are IEA for now.
       my $watch_string = "integral component of";
       my $watch_string_the = "integral component of the";
       my $helper_string_0 = " is localized to the ";
       my $helper_string_1 = " is an ";
       if (@goterms){
         my $count = 0;
         my $verb_count = 0;
         my $goterm="";
        foreach $goterm (@goterms){
#        chomp($goterm);

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
           if ($count gt 0) {
            if ($verb_goterm =~ $helper_string){
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

       } elsif (@goterms gt 2) {
            $helper_string = $helper_string_0;
           if ( $goterm =~/$watch_string/){
                $helper_string = $helper_string_1;
           }
         if ($count eq 0) {
             $verb_goterm = $helper_string . $goterm;
         } elsif ($count lt  @goterms-1) {

            if ($verb_goterm =~ $helper_string){
             $verb_goterm .= "\, the ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= " and "; 
             $verb_goterm .= $helper_string;    
             $verb_goterm .= $goterm;               
            }

         } elsif ($count ge @goterms-1) {
            if ($verb_goterm =~ $helper_string){
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
sub get_gene_protein_xrefs_hash{
 my $species = shift;
 my $xrefs   = shift;
# print "species xrefs $species\t$xrefs\n";
 my %gene_protein_hash=();
 my @lines = read_file($xrefs);
 foreach my $line (@lines){
  my ($gene_id0, $wb_gene_id, $gene_name, $gene_id1, $gene_protein_id, $id0, $id1, $id2) = split(/\t/, $line);
        $gene_protein_id =~ s/^\s+//;
        $gene_protein_id =~ s/\s+$//;
        $wb_gene_id =~ s/^\s+//;
        $wb_gene_id =~ s/\s+$//;
        my $key = $wb_gene_id;
        my $value = $gene_protein_id;
     if ($gene_protein_id ne "\."){
         $gene_protein_hash{$key} = $value;
         }
 }
 return \%gene_protein_hash;
}
sub get_inverse_gene_protein_xrefs_hash{
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
#          print "value is $value\n";
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
#
sub get_gene_name_xrefs_hash{
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
sub get_ortholog_sentence{
 my $gene_id = shift;
 my $prefix  = shift;

 my $sentence="";
 my $gene_name = $gene_name_hash{$gene_id};

 my $ortholog;
 my $ortholog_id;
 my $gene_protein;
 my $elegans_protein_id;

 my @gene_processes = ();
 my @gene_functions = ();
 my @gene_components = ();

 my @GO_processes  = ();
 my @GO_functions  = ();
 my @GO_components = ();
 
 my $function_string = "";
 my $component_string = "";
 my $process_string = "";

 if ($prefix) {

  if ($gene_name =~/$prefix/){
      $ortholog = $gene_name;
      $ortholog =~ s/$prefix\-//gi;
      $ortholog =~ s/^\s+//;
      $ortholog =~ s/\s+$//;
#      print "simple ortholog is $ortholog\n";
      $ortholog_id = $wbgene_elegans_name_hash{$ortholog};
#      print "simple ortholog id is $ortholog_id\n";
  } else {
         $gene_protein = $gene_protein_hash{$gene_id};
     if ($gene_protein){
         $gene_protein =~ s/^\s+//;
         $gene_protein =~ s/\s+$//;
        $elegans_protein_id = $species_ce_protein_hash{$gene_protein};
    }
     if ($elegans_protein_id){
         $elegans_protein_id =~ s/^\s+//;
         $elegans_protein_id =~ s/\s+$//;
        $ortholog = $elegans_name_xrefs_hash{$elegans_protein_id};
        if ($ortholog){
        $ortholog =~ s/^\s+//;
        $ortholog =~ s/\s+$//;
        $ortholog_id = $wbgene_elegans_name_hash{$ortholog};
       }
    }
  }

 }

  if (($ortholog) and ($ortholog_id)){
      $sentence = $ortholog . "\t" . $ortholog_id;
#      print "sentence is $sentence\n"; 

# Check process
 my $gene_process = $elegans_process_hash{$ortholog_id};
 if ($gene_process){
#    print "gene process is $gene_process\n";
  my @ortholog_processes = split(/\,/, $gene_process);
  my $exp_ortholog_processes_ref = get_exp_go_array(\@ortholog_processes);
  my @exp_ortholog_processes = @$exp_ortholog_processes_ref;
     $gene_process = join("\,",@exp_ortholog_processes);
     my $gene_row = $gene_id . "\," . $ortholog . "\," . $gene_process;
     $gene_row =~ s/,+$//;
#     print "gene row is $gene_row\n";
# check for granularity
     my $processes_ref = ConciseDescriptions::get_granular_array($ontology, $gene_row, "IEA", "ISS", \%parents);
     my @processes = @{$processes_ref};

if (@processes){
# convert the ids into terms and put them in the array
    my $process_count = 0;
    foreach my $process (@processes){
     $process_count++;
     if ($process_count gt 2){
        $process =~ s/^\s+//;
        $process =~ s/\s+$//;
# Remove evidence codes
        my ($evidence_to_remove) = $process =~ /\[(.*)\]/;
        if ($evidence_to_remove) {
             $process =~ s/\[[^\]]*\]//g;
           }
     my $go_term = ConciseDescriptions::get_ontology_term($process, \%gene_ontology);    
     push(@GO_processes, $go_term);
     }
    }
     if (@GO_processes){

     my $process_goterm = get_verb_process_goterm(\@GO_processes);
     if ($process_goterm) {   
        $process_string = $process_goterm;
#        print "process_goterm=$process_goterm\n";
        $process_string =~ s/^\s+//;
        $process_string =~ s/\s+$//;
        $process_string =~ s/ \,/\,/g;
        $process_string =~ s/ +/ /g;
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
     }# process_goterm
    }
  }
 }
# Check function
 my $gene_function = $elegans_function_hash{$ortholog_id};
 if ($gene_function){
#  print "gene function is $gene_function\n";
  my @ortholog_functions = split(/\,/, $gene_function);
  my $exp_ortholog_functions_ref = get_exp_go_array(\@ortholog_functions);
  my @exp_ortholog_functions = @$exp_ortholog_functions_ref;
     $gene_function = join("\,",@exp_ortholog_functions);
  my $gene_row = $gene_id . "\," . $ortholog . "\," . $gene_function;
     $gene_row =~ s/,+$//;
#     print "gene row is $gene_row\n";
# check for granularity
  my $functions_ref = ConciseDescriptions::get_granular_array($ontology, $gene_row, "", "", \%parents);
  my @functions = @{$functions_ref};
 if (@functions){
    my $function_count=0;
    foreach my $function (@functions){
          $function_count++;
          if ($function_count gt 2){
# Remove evidence codes
          my ($evidence) = $function =~ /\[(.*)\]/;
             if ($evidence) {
                 $function =~ s/\[$evidence\]//g;
             }
             $function =~ s/^\s+//;
             $function =~ s/\s+$//;
     my $go_term = ConciseDescriptions::get_ontology_term($function, \%gene_ontology);  
     if ( $go_term !~ /activity/) { 
          $go_term =~s/binding/binding activity/g;
     }  
     push(@GO_functions, $go_term);
    }
    } # foreach functions
    if (@GO_functions){
      $function_string = get_verb_function_goterm_exp(\@GO_functions);
#      $function_string .= "\;\n";
      $function_string =~ s/^\s+//;
      $function_string =~ s/\s+$//;
      $function_string =~ s/ \,/\,/g;
      $function_string =~ s/ +/ /g;
      $function_string =~ s/$the_the/$the/g;
      $function_string =~ s/$blank_comma/$comma/g;
      $function_string =~ s/$nadp//g;
      $function_string =~ s/$doublecomma/$comma/g;     
      $function_string =~ s/ +/ /g;
     }
 } # if gene_functions array
}
#
 my $gene_component = $elegans_component_hash{$ortholog_id};
 if ($gene_component){
#  print "gene component is $gene_component\n";
  my @ortholog_components = split(/\,/, $gene_component);
  my $exp_ortholog_components_ref = get_exp_go_array(\@ortholog_components);
  my @exp_ortholog_components = @$exp_ortholog_components_ref;
     $gene_component = join("\,",@exp_ortholog_components);
  my $gene_row = $gene_id . "\," . $ortholog . "\," . $gene_component;
     $gene_row =~ s/,+$//;
#     print "gene row is $gene_row\n";
   my $components_ref = ConciseDescriptions::get_granular_array($ontology, $gene_row, "IEA", "ISS", \%parents);
   my @components = @{$components_ref};
   my $component_count = 0;
   if (@components) {
    $component_count++;
    if ($component_count gt 2){
    foreach my $component (@components){
             my ($evidence) = $component =~ /\[(.*)\]/;
             if ($evidence) {
                 $component =~ s/\[$evidence\]//g;
             }
             $component =~ s/^\s+//;
             $component =~ s/\s+$//;
           
     my $go_term = ConciseDescriptions::get_ontology_term($component, \%gene_ontology);
#     print "go term is $component\t$go_term\n";
     push(@GO_components, $go_term);
    }
   } # foreach components
       if (@GO_components){
        $component_string = get_verb_component_goterm(\@GO_components);

        $component_string =~ s/$component_watch_string/$component_watch_string_the/g;
        $component_string =~ s/$the_the/$the/g;
        $component_string =~ s/$blank_comma/$comma/g;

        $component_string =~ s/ +/ /g;
#        $component_string .= "\n";
        $component_string =~ s/^\s+//;
        $component_string =~ s/\s+$//;
        $component_string =~ s/ \,/\,/g;
        $component_string =~ s/ +/ /g;
      }
   } # if components array
  }
  $sentence = "";
  if ((@GO_components) and (@GO_processes) and (@GO_functions)){
    $sentence = $ortholog . " " . $process_string . " " . $function_string . " and " . $component_string;
  } elsif ((@GO_processes) and (@GO_functions)){
    $sentence = $ortholog . " " . $process_string . " and " . $function_string;
  } elsif ((@GO_components) and (@GO_functions)){
    $sentence = $ortholog . " " . $function_string . " and " . $component_string;
  } elsif ((@GO_components) and (@GO_processes)){
    $sentence = $ortholog . " " . $process_string . " and " . $component_string;
  } elsif (@GO_components){
    $sentence = $ortholog . " " . $component_string;
  } elsif (@GO_processes){
    $sentence = $ortholog . " " . $process_string;
  } elsif (@GO_functions){
    $sentence = $ortholog . " " . $function_string;
  } else { 
    $sentence = "";
  }
# add semi-colon at the end
#    $sentence .= "\;";
#  print "$sentence\n\n";
   if ($sentence){
          my $doublespace = "  ";
          my $space = " ";
# remove any double spaces
          $sentence =~s/$doublespace/$space/g;
     }
# add semi-colon at end
#          $sentence .= "\;";
 }
 return $sentence;
}
sub get_exp_go_array{
 my $go_array_ref = shift;
 my @array=();
 my @go_array = @$go_array_ref;
  foreach my $go (@go_array){
        $go =~ s/^\s+//;
        $go =~ s/\s+$//;
   if ( ($go =~/\[EXP\]/) or 
  ($go =~/\[IDA\]/) or ($go =~/\[IPI\]/) or 
  ($go =~/\[IMP\]/) or ($go =~/\[IMP\]/) or 
  ($go =~/\[IGI\]/) or ($go =~/\[IEP\]/) ) {
    push(@array, $go);  
  }
 }
  my @sorted_array = sort(uniq(@array));
  
 return \@sorted_array;
}
