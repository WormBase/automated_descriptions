#!/usr/bin/env perl
use ConciseDescriptions;
use LWP::Simple;
use List::MoreUtils qw(uniq);
use List::Util qw/first/;
use File::Slurp;
use Text::CSV;
use warnings;
use strict;

my $AND = "AND";
my $species_term = $ARGV[0];
my ($species, $project, $species_name, $species_prefix) = split(/$AND/, $species_term);

chomp($species);
chomp($project);
chomp($species_name);
chomp($species_prefix);
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

my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
my $gene_list_dir = $home . "release/$PRODUCTION_RELEASE/c_elegans/gene_lists/";
my $db_gene_list  = $gene_list_dir . "wb_gene_list.txt";
#
# Define component adjectives
my @component_adjectives = ("GO\:0005622");
#
# Define GO term file
my $ontology = "GO";
my $go_file = $gene_list_dir . "go_terms.txt";
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
my $gene_association_elegans_file = $go_elegans_dir . "gene_association.wb";
my $orthology = $home . "release/$PRODUCTION_RELEASE/$species/orthology/";
my $elegans_orthology = $home . "release/$PRODUCTION_RELEASE/c_elegans/orthology/";

my $go_path = $home . "release/$PRODUCTION_RELEASE/$species/gene_ontology/output_files/";
my $output_file = $go_path . "sentences_for_GO.txt";
#
my $doublecomma = "\,\,";
my $the_the = "the the";
my $the = "the";
my $is_is = "is is";
my $is = "is";
my $is_a_is_a = "is a is a";
my $is_a = "is a";
my $blank_comma = "\, \,";
my $comma = "\, ";
my $nadp =" NAD or NADP as acceptor\,";
my $component_watch_string = "integral component of";
my $component_watch_string_the = "integral component of the";
#my $is_localized_cellular="is localized to the intracellular";
#my $is_cellular="is intracellular";
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

my $the_cell_0 = "localized to the cell";
my $the_cell_1 = "expressed widely";

my $growth_0 = "multicellular organism growth";
my $growth_1 = "growth";
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
my $gene_name_hash_ref = get_gene_name_hash($gene_association_species_file, $db_gene_list);
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
foreach my $test (@live_gene_array){
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
 print "$gene_id\n";
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
     next if ($process=~/$gene_id/);
     next if ($process=~/$gene_name/);

     if ($process_count > 0){
        $process =~ s/^\s+//;
        $process =~ s/\s+$//;
     my $go_term = ConciseDescriptions::get_ontology_term($process, \%gene_ontology); 
     if ($go_term){
      if (lc($go_term) !~/obsolete/){   
        push(@GO_processes, $go_term);
       }
      }
#     print "process go_term $process\t$go_term\n";
     }
    }
#
     if (@GO_processes){
     @GO_processes=uniq(sort(@GO_processes));
     my $process_goterm = get_verb_process_goterm(\@GO_processes);
     if ($process_goterm) {   
        $process_string = $process_goterm;
        print "process $process_string\n";
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
 print "function\: $gene_row\n";
# check for granularity
 my $functions_ref = ConciseDescriptions::get_granular_array($ontology, $gene_row, "", "", \%parents);
 my @functions = @{$functions_ref};
 if (@functions){
    my $function_count=0;
    foreach my $function (@functions){
          $function_count++;
     next if ($function=~/$gene_id/);
     next if ($function=~/$gene_name/);

          if ($function_count > 0){
          my ($evidence) = $function =~ /\[(.*)\]/;
             if ($evidence) {
                 $function =~ s/\[$evidence\]//g;
             }
             $function =~ s/^\s+//;
             $function =~ s/\s+$//;
     my $go_term = ConciseDescriptions::get_ontology_term($function, \%gene_ontology);  
     if ($go_term){
     if ( $go_term !~ /activity/) { 
          $go_term =~s/binding/binding activity/g;
     }
     if ( lc($go_term) !~/obsolete/ ) {
       push(@GO_functions, $go_term);
     }
    }
    }
    } # foreach functions
    if (@GO_functions){
      @GO_functions = uniq(sort(@GO_functions));
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
    @gene_components = split(/\,/, $gene_component);
 }
   my $components_ref = ConciseDescriptions::get_granular_array($ontology, $gene_row, "", "", \%parents);
   my @components = @{$components_ref};
#
# Check if any components are adjectives 
#
   my @new_component_adjectives=();
   if (@components){
     my $pointer = 0;
      while ($pointer <= $#components) {
       my $component = $components[$pointer];
       foreach my $adjective (@component_adjectives){
        if ($component =~/$adjective/){
           push(@new_component_adjectives, $component);
           splice(@components, $pointer, 1);
        } else {
           $pointer++;
        } 
       } 
     }
   }

 my $adj_phrase = " is ";
 my $adj_sentence;
 my $adj_count = @new_component_adjectives; 

 if ($adj_count > 0){
   my $count = 0;
   $adj_sentence = $adj_phrase;
   foreach my $adj (@new_component_adjectives){
             my ($evidence) = $adj =~ /\[(.*)\]/;
             if ($evidence) {
                 $adj =~ s/\[$evidence\]//g;
             }
             $adj =~ s/^\s+//;
             $adj =~ s/\s+$//;
     my $go_term = "";
     $go_term = ConciseDescriptions::get_ontology_term($adj, \%gene_ontology);
     next if ($go_term eq "");
     next if (lc($go_term) =~/obsolete/);
     my $count++;
     if (($count == 1) and ($adj_count < 3)) {
      $adj_sentence .= $go_term;
     } elsif (($count < ($adj_count-1)) and ($adj_count > 2)){
       $adj_sentence .= $go_term . "\,";
     } elsif (($count == $adj_count) and ($adj_count > 1)){
       $adj_sentence .= " and " . $go_term;
     }
  }
 }

#
   my $component_count = 0;
   if (@components) {
    foreach my $component (@components){
    $component_count++;
     next if ($component=~/$gene_id/);
     next if ($component=~/$gene_name/);

    if ($component_count > 0){
             my ($evidence) = $component =~ /\[(.*)\]/;
             if ($evidence) {
                 $component =~ s/\[$evidence\]//g;
             }
             $component =~ s/^\s+//;
             $component =~ s/\s+$//;
     my $go_term = "";
        $go_term = ConciseDescriptions::get_ontology_term($component, \%gene_ontology);
     if ($go_term){
      if (lc($go_term) !~/obsolete/) {  
       push(@GO_components, $go_term);
      }
     }
    }
   } # foreach components
       if (@GO_components){
        @GO_components = uniq(sort(@GO_components));
        $component_string = get_verb_component_goterm(\@GO_components);
        print "component $component_string\n";
        $component_string =~ s/$component_watch_string/$component_watch_string_the/g;
        $component_string =~ s/$structural/$is_structural/g;
        $component_string =~ s/$the_the/$the/g;
        $component_string =~ s/$blank_comma/$comma/g;

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
          $sentence =~ s/$structural/$is_structural/g;
          $sentence =~ s/$is_is/$is/g;
          $sentence =~s/$is_a_is_a/$is_a/g;
          $sentence =~s/$doublespace/$space/g;
          $sentence =~s/$the_cell_0/$the_cell_1/g;
          $sentence =~s/$growth_0/$growth_1/g;
# if sentence has data based on protein domain information, that information should be first.
         if ($species !~/elegans/){
          if ($sentence =~/$helper_string_after_iea/){
           my $tmp = $sentence;
           $tmp =~s/$helper_string_after_iea//g;
           $sentence = $helper_string_before_iea . $tmp;
          }
          if (($sentence !~/$helper_string_before_iea/) and ($sentence)) {
           my $tmp = $sentence;
           $sentence = $helper_string_before_iea . $tmp;
          }
         }
# add semi-colon at end
#          $sentence .= "\;";
# add ortholog sentence
           my $ortholog_sentence = "";
      if (length($sentence) > 0){
              $sentence .= "\.";
      
       my $out = $individual_path . $gene_id;
        write_file($out, $sentence);
        write_file($output_file, {append => 1 }, $sentence);
        write_file($output_file, {append => 1 }, "\n\n\n");
       }
 }
exit(0);

sub get_gene_array{
 my $file = shift;
 my @array=();
 my @lines = read_file($file);
  foreach my $line (@lines){
   next if ($line =~/\!/);
   chomp($line);
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
 my $second_file = shift;
 my %hash=();
 my @lines = read_file($file);
  foreach my $line (@lines){
   next if ($line =~/\!/);
   chomp($line);
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
 my @second_lines = read_file($second_file);
  foreach my $line (@second_lines){
   chomp($line);
   my @fields = split(/\t/, $line);
   my $gene_id = $fields[0];
   my $gene_name = $fields[1];
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
   my $GO = $fields[8]; 
   if ($ec){
       $term = $go_term . "\[" . $evidence_code . "\]";
   } else{
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
         my ($evidence) = $temp_string =~ /$go_term\[(\w+)\]/;
         print "evidence is $evidence\n";
         if (grep {$_ =~ $evidence} @exp_ec) {
          print "$gene_id\t$go_term\t$evidence\t$hash{$gene_id}\n";
         } elsif ($evidence !~ $evidence_code) {
          my $old =  $go_term . "\[" . $evidence . "\]";
          my $replace = $go_term . "\[" . $evidence_code . "\]";
            $temp_string =~ s/$old/$replace/g;
            $hash{$gene_id} = $temp_string;
         }
#         my $replace = $go_term . "\[" . $evidence_code . "\]";
#            $temp_string =~ s/$go_term\[[^\]]*\]/$replace/g;
#            $hash{$gene_id} = $temp_string;
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
           if ($count > 0) {
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

       } elsif (@goterms > 2) {
            $helper_string = $helper_string_0;
           if ( $goterm =~/$watch_string/){
                $helper_string = $helper_string_1;
           }
         if ($count eq 0) {
             $verb_goterm = $helper_string . $goterm;
         } elsif ($count < (@goterms-1)) {

            if ($verb_goterm =~ $helper_string){
             $verb_goterm .= "\, ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= " and "; 
             $verb_goterm .= $helper_string;    
             $verb_goterm .= $goterm;               
            }

         } elsif ($count >= (@goterms-1)) {
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

           if ($count > 0) {

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
         } elsif ($count >= (@goterms-1)) {
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

           if ($count > 0) {

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
         } elsif ($count >= (@goterms-1)) {
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
           if ($count > 0) {
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

       } elsif (@goterms > 2) {
            $helper_string = $helper_string_0;
           if ( $goterm =~/$watch_string/){
                $helper_string = $helper_string_1;
           }
         if ($count == 0) {
             $verb_goterm = $helper_string . $goterm;
         } elsif ($count < (@goterms-1)) {

            if ($verb_goterm =~ $helper_string){
             $verb_goterm .= "\, the ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= " and "; 
             $verb_goterm .= $helper_string;    
             $verb_goterm .= $goterm;               
            }

         } elsif ($count >= (@goterms-1)) {
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
