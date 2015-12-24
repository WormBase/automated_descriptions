#!/usr/bin/env perl
#
use strict;
use warnings;
use List::MoreUtils qw(uniq);
use Switch;
use Spreadsheet::XLSX;
use File::Slurp;
use File::Copy;
use Data::Dumper;
use ConciseDescriptions;
#
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
my $html = ConciseDescriptions::get_html_dir();
my $home_elegans = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/c_elegans/";
my $tissue_dir = $home_elegans . "tissue_expression/";
my $gene_list_dir = $home_elegans . "gene_lists/" ;
my $input_path = $tissue_dir . "input_files/";
my $old_input_path = $html . "concise_descriptions/release/$RELEASE/c_elegans/tissue_expression/input_files/";
my $xlsx_file = $input_path . "Neuron_list_Ranjana.xlsx";
my $old_xlsx_file = $old_input_path . "Neuron_list_Ranjana.xlsx";
copy($old_xlsx_file, $xlsx_file);
#
my $vnc = "ventral cord neuron"; 
my $vnc2 = "vnc motor neuron"; 
my $interneuron = "interneuron";
my $sensory = "sensory neuron";
my $head = "head motor neuron";
my $motor = "motor neuron";
my $mechanosensory = "mechanosensory neuron";
my $inner_labial = "inner labial neuron";
my $outer_labial = "outer labial neuron";
my $pn = "pharyngeal neuron";
my $pmn = "pharyngeal motorneuron";
my $pmns = "pharyngeal motorneurons";
my $amphid = "amphid sensory neuron";
my $pharyngeal_interneuron = "pharyngeal interneuron";

my $output_dir = $input_path;
my $vnc_output = $output_dir . "ventral_cord_neuron.txt";
my $interneuron_output = $output_dir . "interneuron.txt";
my $sensory_output = $output_dir . "sensory_neuron.txt";
my $head_output =  $output_dir . "head_neuron.txt";
my $motor_output = $output_dir . "motor_neuron.txt";
my $mechanosensory_output = $output_dir . "mechanosensory_neuron.txt";
my $inner_labial_output = $output_dir . "inner_labial_neuron.txt";
my $outer_labial_output = $output_dir . "outer_labial_neuron.txt";
my $pn_output = $output_dir . "pharyngeal_neuron.txt";
my $pmn_output = $output_dir . "pharyngeal_motor_neuron.txt";
my $amphid_output = $output_dir . "amphid.txt";
my $pharyngeal_interneuron_output = $output_dir . "pharyngeal_interneuron.txt";

if (-e $vnc_output){
   my @args =("/bin/rm", "-f", $vnc_output);
   system(@args) == 0 or die("could not delete file $vnc_output\n");
}
if (-e $interneuron_output){
   my @args =("/bin/rm", "-f", $interneuron_output);
   system(@args) == 0 or die("could not delete file $interneuron_output\n");
}
if (-e $sensory_output){
   my @args =("/bin/rm", "-f", $sensory_output);
   system(@args) == 0 or die("could not delete file $sensory_output\n");
}
if (-e $head_output){
   my @args =("/bin/rm", "-f", $head_output);
   system(@args) == 0 or die("could not delete file $head_output\n");
}
if (-e $motor_output){
   my @args =("/bin/rm", "-f", $motor_output);
   system(@args) == 0 or die("could not delete file $motor_output\n");
}
if (-e $mechanosensory_output){
   my @args =("/bin/rm", "-f", $mechanosensory_output);
   system(@args) == 0 or die("could not delete file $mechanosensory_output\n");
}
if (-e $inner_labial_output){
   my @args =("/bin/rm", "-f", $inner_labial_output);
   system(@args) == 0 or die("could not delete file $inner_labial_output\n");
}
if (-e $outer_labial_output){
   my @args =("/bin/rm", "-f", $outer_labial_output);
   system(@args) == 0 or die("could not delete file $outer_labial_output\n");
}
if (-e $pmn_output){
   my @args =("/bin/rm", "-f", $pmn_output);
   system(@args) == 0 or die("could not delete file $pmn_output\n");
}
if (-e $pn_output){
   my @args =("/bin/rm", "-f", $pn_output);
   system(@args) == 0 or die("could not delete file $pn_output\n");
}
if (-e $amphid_output){
   my @args =("/bin/rm", "-f", $amphid_output);
   system(@args) == 0 or die("could not delete file $amphid_output\n");
}
if (-e $pharyngeal_interneuron_output){
   my @args =("/bin/rm", "-f", $pharyngeal_interneuron_output);
   system(@args) == 0 or die("could not delete file $pharyngeal_interneuron_output\n");
}
my @not_in_ontology_array=();
my @unique_not_in_ontology_array=();
# Create WBbt hash so that terms can be referenced by WBbt ID.
my %anatomy_ontology;
my $anatomy_file = $input_path . "anatomy_terms.txt";
my @anat_lines = read_file($anatomy_file); 
my $synonym_file = $input_path . "anatomy_synonyms.txt";
my @syn_lines = read_file($synonym_file);
foreach my $syn_line (@syn_lines){
 chomp($syn_line);
 my $key   = $syn_line;
 my $value = $syn_line;
 $value  =~ /\((WBbt.+?)\)/;
 $value = $1; 
 $key =~ s/\(.*//g;
 $key =~ s/^\s+//;
 $key =~ s/\s+$//;
 chomp($key);
 chomp($value);
 $anatomy_ontology{$key} = $value;
}
foreach my $anat_line (@anat_lines){
 chomp($anat_line);
 my $key   = $anat_line;
 my $value = $anat_line;
 $value  =~ /\((WBbt.+?)\)/;
 $value = $1; 
 $key =~ s/\(.*//g;
 $key =~ s/^\s+//;
 $key =~ s/\s+$//;
 chomp($key);
 chomp($value);
 if ($value =~/WBbt\:0005339/){
  print "WBbt\:0005339 $key\t$value\n"
 }
 $anatomy_ontology{$key} = $value;
}
#
my $excel = Spreadsheet::XLSX -> new ($xlsx_file);
 
foreach my $sheet (@{$excel -> {Worksheet}}) {
 
       printf("Sheet: %s\n", $sheet->{Name});
        
       $sheet -> {MaxRow} ||= $sheet -> {MinRow};
        
        foreach my $row ($sheet -> {MinRow} .. $sheet -> {MaxRow}) {
         
               $sheet -> {MaxCol} ||= $sheet -> {MinCol};
                my $zero_column = "";
                my $first_column = "";
                my $second_column = "";
                my $third_column = "";
                my $ROW = $row+1;
                my $category = "";
                my $COL="";
               foreach my $col ($sheet -> {MinCol} ..  $sheet -> {MaxCol}) {
                       my $COL = $col;
                       my $cell = $sheet -> {Cells} [$row] [$col];
                       my $zero_column_key = $sheet -> {Cells} [$row] [0] -> {Val};
                       if ($zero_column_key){
                        chomp($zero_column_key);
                        $zero_column_key =~ s/^\s+//;
                        $zero_column_key =~ s/\s+$//;
                       }
                       if ($zero_column_key) {
                           print "zero column key $zero_column_key\n";                           
                         if ($anatomy_ontology{$zero_column_key}) {
                           $zero_column = "$zero_column_key\t$anatomy_ontology{$zero_column_key}\n";
                           print "zero column $zero_column";
                        } else{
                           push(@not_in_ontology_array, $zero_column_key);
                        }
                       }
                       my $first_column_key = $sheet -> {Cells} [$row] [1] -> {Val};
                       if ($anatomy_ontology{$first_column_key}) {
                        $first_column = "$first_column_key\t$anatomy_ontology{$first_column_key}\n";
                       } else {
                          push(@not_in_ontology_array, $first_column_key);
                       }
                       my $second_column_key = $sheet -> {Cells} [$row] [2] -> {Val};
                       if ($anatomy_ontology{$second_column_key}) {
                         $second_column = "$second_column_key\t$anatomy_ontology{$second_column_key}\n";
                       } else {
                          push(@not_in_ontology_array, $second_column_key);
                       }
                       $third_column = $sheet -> {Cells} [$row] [3] -> {Val};
                       $category = lc($third_column);

                       if ($cell) {
                            my $category = "";
                            my $value = lc $cell -> {Val};
                            if ($col==3) {
                                $category = $value;
                           }
                           
                       }
 
               } # column
     			switch ($category) {
        			case ($vnc) {  write_output($vnc_output, $zero_column, $first_column, $second_column); }
        			case ($vnc2) { write_output($vnc_output, $zero_column, $first_column, $second_column); }
        			case ($interneuron) { 
                                               write_output($interneuron_output, $zero_column, $first_column, $second_column); }
        			case ($sensory) { 
                                               write_output($sensory_output, $zero_column, $first_column, $second_column); }
        			case ($mechanosensory) { 
                                               write_output($mechanosensory_output, $zero_column, $first_column, $second_column); }
        			case ($head) { 
                                               write_output($head_output, $zero_column, $first_column, $second_column); }
        			case ($motor) { 
                                               write_output($motor_output, $zero_column, $first_column, $second_column); }
        			case ($amphid) { 
                                               write_output($amphid_output, $zero_column, $first_column, $second_column); }
        			case ($pn) { 
                                               write_output($pn_output, $zero_column, $first_column, $second_column); }
        			case ($pmn) { 
                                               write_output($pmn_output, $zero_column, $first_column, $second_column); }
        			case ($pmns) { 
                                               write_output($pmn_output, $zero_column, $first_column, $second_column); }
        			case ($inner_labial) { 
                                               write_output($inner_labial_output, $zero_column, $first_column, $second_column); }
        			case ($outer_labial) { write_output($outer_labial_output, $zero_column, $first_column, $second_column);                                                   }
        			case ($pharyngeal_interneuron) { 
                                              write_output($pharyngeal_interneuron_output, $zero_column, $first_column, $second_column); }
		        else  { print "no output file for $category\t row\:$ROW\n"; }
			    }
       } # row
 
}

@unique_not_in_ontology_array = sort(uniq(@not_in_ontology_array));
foreach my $element (@unique_not_in_ontology_array){
 print "$element not in ontology\.\n";
}

exit 0;

sub write_output {
 my $output = shift;
 my $zero_column = shift;
 my $first_column = shift;
 my $second_column = shift;

   if ($zero_column) {
        write_file($output, {append => 1}, $zero_column);
       }
   if ($first_column) {
        write_file($output, {append => 1}, $first_column);
       }
   if ($second_column) {  
        write_file($output, {append => 1}, $second_column);
     } 
}
