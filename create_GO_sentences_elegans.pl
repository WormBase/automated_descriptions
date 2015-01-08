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

if ($species =~/elegans/){
#
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/c_elegans/";
my $gene_list_dir = $home . "gene_lists/";
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
my $go_elegans_dir = $home . "gene_ontology/input_files/";
my $gene_association_elegans_file = $go_elegans_dir . "gene_association.$RELEASE.wb.c_elegans";
my $elegans_orthology = $home . "orthology/";
my $go_path = $home . "gene_ontology/output_files/";
my $output_file = $go_path . "sentences_for_GO.txt";
#
my $wbgene_elegans_name_hash_ref = get_wbgene_name_elegans_hash($gene_association_elegans_file);
my %wbgene_elegans_name_hash = %$wbgene_elegans_name_hash_ref;
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
if (-e $go_path){
   my $go_path_flag = 1;
#   print "$go_path exists\n"; 
} else {
   mkdir $go_path or die "could not create $go_path";
}
my $individual_path = $go_path . "individual_gene_sentences/";
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

#my $infile = $gene_list_dir . "sort.wbgene_list.no_concise_descriptions.txt";
my $infile = $gene_list_dir . "sort.uncurated_genes.txt";

my @gene_array = read_file($infile);
my $gene_name_hash_ref = get_gene_name_hash($gene_association_elegans_file);
my %gene_name_hash = %$gene_name_hash_ref;
my $go = "P";
my $elegans_process_hash_ref = get_gene_go_hash($gene_association_elegans_file, $go, 1);
my %elegans_process_hash = %$elegans_process_hash_ref;
   $go = "C";
my $elegans_component_hash_ref = get_gene_go_hash($gene_association_elegans_file, $go, 1);
my %elegans_component_hash = %$elegans_component_hash_ref;
   $go = "F";
my $elegans_function_hash_ref = get_gene_go_hash($gene_association_elegans_file, $go, 1);
my %elegans_function_hash = %$elegans_function_hash_ref;

foreach my $gene_id (@gene_array){
 chomp($gene_id);
       $gene_id =~ s/^\s+//;
       $gene_id =~ s/\s+$//;

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
 my $gene_process = $elegans_process_hash{$gene_id};
 if ($gene_process){
 my $gene_row = $gene_id . "\," . $gene_name . "\," . $gene_process;
# print "gene row is $gene_row\n";
# check for granularity
 my $processes_ref = ConciseDescriptions::get_granular_array($ontology, $gene_row, "IEA", "ISS", \%parents);
 my @processes = @{$processes_ref};
 if (@processes){
# convert the ids into terms and put them in the array
    my $process_count = 0;
    foreach my $process (@processes){
     next if ($process =~/\[IEA\]/);
     next if ($process =~/\[ISS\]/);
     $process_count++;
     if ($process_count gt 2){
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
         $go_term .= "[$evidence]";
         push(@GO_processes, $go_term);
      }
     }
    }
# 
#     my $process_goterm = join(', ', @GO_processes)
#
     if (@GO_processes){
#
# yet untested, $imp="IMP" would identify processes with evidence code separately; else it includes them without 
# adding information that it was "based on mutant phenotype"
#
     my $imp = "";
     my $process_goterm = get_verb_process_goterm(\@GO_processes, $imp);
     if ($process_goterm) {   
        $process_string = $process_goterm;
#        print "process string = $process_string\n";
#        $process_string .= "\;\n";
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
     }# process_goterm
    }
 } # gene_process
}
 my $gene_function = $elegans_function_hash{$gene_id};
 if ($gene_function){
 my $gene_row = $gene_id . "\," . $gene_name . "\," . $gene_function;
# check for granularity
 my $functions_ref = ConciseDescriptions::get_granular_array($ontology, $gene_row, "", "", \%parents);
 my @functions = @{$functions_ref};
 my $exp_functions_ref = get_exp_go_array(\@functions);
 my @exp_functions=@$exp_functions_ref;
 my $iss_functions_ref = get_stat_go_array(\@functions);
 my @iss_functions=@$iss_functions_ref;
 my $iea_functions_ref = get_iea_go_array(\@functions);
 my @iea_functions=@$iea_functions_ref;
 my $iea_function_string="";
 my $iss_function_string="";
 my $exp_function_string="";
 
 if (@iea_functions){
      my $ieaf_string="";
      foreach my $ieaf (@iea_functions){
       $ieaf_string .= "$ieaf\t";
     }
     $iea_function_string = get_function_string($iea_functions_ref, "IEA", \%gene_ontology);
#     print "iea string\t$ieaf_string\t$iea_function_string\n";
 }
 if (@iss_functions){
      my $issf_string="";
      foreach my $issf (@iss_functions){
       $issf_string .= "$issf\t";
     }
     $iss_function_string = get_function_string($iss_functions_ref, "ISS", \%gene_ontology);
#     print "iss string\t$issf_string\t$iss_function_string\n";
 }
 if (@exp_functions){
     $exp_function_string = get_function_string($exp_functions_ref, "EXP", \%gene_ontology);
 }
 foreach my $function (@iea_functions){
  push (@GO_functions, $function);
 }
 foreach my $function (@iss_functions){
  push (@GO_functions, $function);
 }
 foreach my $function (@exp_functions){
  push (@GO_functions, $function);
 }
 if ((@exp_functions) and (@iss_functions) and (@iea_functions)){
     $function_string = $exp_function_string . " and " . $iss_function_string . " and " . $iea_function_string;
 } elsif ((@exp_functions) and (@iss_functions)) {
     $function_string = $exp_function_string . " and " . $iss_function_string;
 } elsif ((@exp_functions) and (@iea_functions)) {
     $function_string = $exp_function_string . " and " . $iea_function_string;
 } elsif ((@iss_functions) and (@iea_functions)) {
     $function_string = $iss_function_string . " and " . $iea_function_string;
 } elsif (@exp_functions) {
     $function_string = $exp_function_string;
 } elsif (@iss_functions) {
     $function_string = $iss_function_string;
 } elsif (@iea_functions){
     $function_string = $iea_function_string;
 }
#      my ($evidence) = $function_string =~ /\[(.*)\]/;
#      if ($evidence) {
#          $function_string =~ s/\[$evidence\]//g;
#      }
      $function_string =~ s/^\s+//;
      $function_string =~ s/\s+$//;
      $function_string =~ s/ \,/\,/g;
      $function_string =~ s/ +/ /g;
      $function_string =~ s/$the_the/$the/g;
      $function_string =~ s/$blank_comma/$comma/g;
      $function_string =~ s/$nadp//g;
      $function_string =~ s/$doublecomma/$comma/g;     
      $function_string =~ s/ +/ /g;

} # if gene function string
 my $gene_component = $elegans_component_hash{$gene_id};
 if ($gene_component){
  my $gene_row = $gene_id . "\," . $gene_name . "\," . $gene_component;
 if ($gene_component){
    @gene_components = split(/\t/, $gene_component);
 }
   my $components_ref = ConciseDescriptions::get_granular_array($ontology, $gene_row, "IEA", "ISS", \%parents);
   my @components = @{$components_ref};
   my $component_count = 0;
   if (@components) {
    foreach my $component (@components){
    $component_count++;
    if ($component_count gt 2){
     next if ($component =~/\[IEA\]/);
     next if ($component =~/\[ISS\]/);
             my ($evidence) = $component =~ /\[(.*)\]/;
             if ($evidence) {
                 $component =~ s/\[$evidence\]//g;
             }
             $component =~ s/^\s+//;
             $component =~ s/\s+$//;
           
     my $go_term = ConciseDescriptions::get_ontology_term($component, \%gene_ontology);
        $go_term .="[$evidence]";  
     push(@GO_components, $go_term);
    }
   } # foreach components
       if (@GO_components){
        $component_string = get_verb_component_goterm(\@GO_components);

        $component_string =~ s/$component_watch_string/$component_watch_string_the/g;
        $component_string =~ s/$the_the/$the/g;
        $component_string =~ s/$blank_comma/$comma/g;
my @evidences=("EXP","IMP","IGI","IPI","IDA","IEP","IEA","ISS","ISA","ISO","ISM","IGC","IBA","IBD","IKR","IRD","RCA","TAS","NAS","ND","IC");
# Remove evidence codes
         foreach my $ec (@evidences){
	    $component_string =~ s/\[$ec\]//g;
         }
        $component_string =~ s/ +/ /g;
        $component_string =~ s/^\s+//;
        $component_string =~ s/\s+$//;
        $component_string =~ s/ \,/\,/g;
        $component_string =~ s/ +/ /g;
      }
   } # if components array
 }
  my $sentence = "";
  if ((@GO_components) and (@GO_processes) and (@GO_functions)){
    $sentence = $gene_name . " " . $process_string . "\; " . $gene_name . " " . $function_string . "\; " . $gene_name . " " . $component_string;
  } elsif ((@GO_processes) and (@GO_functions)){
    $sentence = $gene_name . " " . $process_string . "\; " . $gene_name . " " . $function_string;
  } elsif ((@GO_components) and (@GO_functions)){
    $sentence = $gene_name . " " . $function_string . "\; " . $gene_name . " " . $component_string;
  } elsif ((@GO_components) and (@GO_processes)){
    $sentence = $gene_name . " " . $process_string . "\; " . $gene_name . " " . $component_string;
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
# add semi-colon at end
#          $sentence .= "\;";
      if (length($sentence) gt 0){
              $sentence .= "\.";
       my $out = $individual_path . $gene_id;
        write_file($out, $sentence);
        write_file($output_file, {append => 1 }, $sentence);
        write_file($output_file, {append => 1 }, "\n\n\n");
       }
 }
}
exit(0);

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
    my $imp = shift;
    my $imp_string = " based on mutant phenotype";
    my @goterms = @$goterms_ref;
    my $verb_goterm="";
# 
# The following was added (untested) to incorporate an IMP, "based on mutant phenotype".
# cluster the GO terms

   if ($imp !~/IMP/){
       $verb_goterm = get_verb_goterm_process(\@goterms);
 } else {     
    my @IMP=();
    my @other_ec =();
    if (@goterms){
       foreach my $goterm (@goterms){
         if ($goterm =~/\[IMP\]/){
            push(@IMP, $goterm);
         } else {
            push(@other_ec, $goterm);
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
    my $helper_string_0 = " is involved in ";
    my $helper_string_1 = " functions in ";
    my $watch_string = "involved in";
    my @goterms = @$goterms_ref;
    my $helper_string="";
    my $verb_goterm="";
#
       if (@goterms){
         my $count = 0;
         my $verb_count = 0;
        foreach my $goterm (@goterms){
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
# print "calculating $verb_goterm\n";
 return $verb_goterm;
}
sub get_verb_component_goterm{
       my $goterms_ref = shift;
       my @goterms = @$goterms_ref;
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
sub get_stat_go_array{
 my $go_array_ref = shift;
 my @array=();
 my @go_array = @$go_array_ref;
  foreach my $go (@go_array){
        $go =~ s/^\s+//;
        $go =~ s/\s+$//;
   if ( ($go =~/\[ISS\]/) or 
  ($go =~/\[ISO\]/) or ($go =~/\[ISA\]/) or 
  ($go =~/\[ISM\]/) or ($go =~/\[IGC\]/) or 
  ($go =~/\[IBA\]/) or ($go =~/\[IBD\]/) or 
  ($go =~/\[IKR\]/) or ($go =~/\[IRD\]/) or
  ($go =~/\[RCA\]/) ) {
    push(@array, $go);  
  }
 }
  my @sorted_array = sort(uniq(@array));
  
 return \@sorted_array;
}
sub get_iea_go_array{
 my $go_array_ref = shift;
 my @array=();
 my @go_array = @$go_array_ref;
  foreach my $go (@go_array){
        $go =~ s/^\s+//;
        $go =~ s/\s+$//;
   if (($go =~/\[IEA\]/) or 
  ($go =~/\[TAS\]/) or ($go =~/\[NAS\]/) or 
  ($go =~/\[IC\]/) or ($go =~/\[ND\]/)) {
    push(@array, $go);  
  }
 }
  my @sorted_array = sort(uniq(@array));
  
 return \@sorted_array;
}
sub get_function_string {
 my $functions_array_ref = shift;
 my $ec = shift;
 my $gene_ontology_ref = shift;
 my $function_string = "";
 my @functions = @$functions_array_ref;
 my @GO_functions = ();
 
 if (@functions){
    my $function_count=0;
    foreach my $function (@functions){
#          print "function is $function\n";
          next if ($function =~/GO\:0005488/); # ignore binding
          next if ($function =~/GO\:0005515/); # ignore protein binding
          $function_count++;
          if ($function_count gt 0){
#          print "function is $function\n";
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
#     print "goterm is $go_term\n";
     push(@GO_functions, $go_term);
    }
    } # foreach functions

    if (@GO_functions){
       my $gof_string = "";
       foreach my $gof (@GO_functions){
        $gof_string .= "$gof\t";
       }
       $function_string = get_verb_function_goterm(\@GO_functions, $ec);
#       print "$gof_string\t\nfunction string is $function_string\n";
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
# print "$size\tfunction string is $function_string\n";

 } # if gene_functions array

 return $function_string;
}
