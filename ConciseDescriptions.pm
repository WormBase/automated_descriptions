package ConciseDescriptions;
use LWP::UserAgent;
use File::Slurp; 
use List::MoreUtils qw(uniq);
use List::Util qw/first/;
use List::MoreUtils qw( minmax );
use Switch;
use Text::CSV;
use strict;
use warnings;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(get_production_release get_release get_html_dir get_cgi_dir get_db_ip getwebpage get_granular_array get_parents get_ontology_term get_array get_ontology_hash get_ontology_parents_children get_species_ensp_hash get_gene_protein_hash);
#
# This script is used to support automated concise descriptions for wormbase:
# http://wiki.wormbase.org/index.php/Generation_of_automated_descriptions 
# J. Done and R. Kishore, California Institute of Technology, 2014.
#
#
sub get_release{
# get_release returns the release of the input data
 my $html = ConciseDescriptions::get_html_dir();
 my $home = $html . "concise_descriptions/";
 my $release_file = $home . "release.txt";
 my $release = read_file($release_file);
 chomp($release);
 $release =~s/^\s+//;
 $release =~s/\s+$//;
 my $RELEASE = uc($release);
 return $RELEASE;
}
sub get_production_release{
# get_production_release returns the production release to which 
# the output of these routines are contributing.
 my $html = ConciseDescriptions::get_html_dir();
 my $home = $html . "concise_descriptions/";
 my $release_file = $home . "production_release.txt";
 my $release = read_file($release_file);
 chomp($release);
 $release =~s/^\s+//;
 $release =~s/\s+$//;
 my $RELEASE = uc($release);
 return $RELEASE;
}
# html path is defined by reading the an input file (e.g., docroot in RH, /var/www in Ubuntu, etc.)
sub get_html_dir {
 use File::Slurp;
 my $html_file = "./html.txt";
 my $html = read_file($html_file);
 chomp($html);
 $html =~s/\s+$//;
 $html =~s/^\s+//;

 return $html;
}
# cgi path is defined by reading the an input file (e.g., cgi-bin in RH, etc.)
sub get_cgi_dir {
 use File::Slurp;
 my $cgi_file = "./cgi.txt";
 my $cgi = read_file($cgi_file);
 chomp($cgi);
 $cgi =~s/\s+$//;
 $cgi =~s/^\s+//;

 return $cgi;
}
sub get_db_ip {
 use File::Slurp;
 my $ip_file = "./db_ip.txt";
 my $dbip;
 $dbip = read_file($ip_file);
 chomp($dbip);
 $dbip =~ s/^\s+//;
 $dbip =~ s/\s+$//;
 return $dbip;
}
sub getwebpage{
    my $u = shift;
    my $page = "";
    use LWP::UserAgent;
	
    my $ua = LWP::UserAgent->new(timeout => 60); # instantiates a new user agent
    my $request = HTTP::Request->new(GET => $u); # grabs url
    my $response = $ua->request($request);       # checks url, dies if not valid.
    if ($response->is_success) {
        $page = $response->content;    #splits by line
    } else {
        warn $response->status_line,"\n";
        $page = $response->status_line;
    }

    return $page;
}
# get_granular_array returns an array without ontology parents leaving the most granular terms
#                    inputs are $ontology_code (e.g., "GO")
#                               $row from processed gene association file.
#                               $evidence_exclusion_1 and $evidence_exclusion_2 (optional) 
#                                                     represent ignored evidence codes
#                               $parents_ref is the reference of the parents hash for each ontology term 
#                     The output is the reference to an array containing the ontology identification numbers 
#                                without the ontology parents.
#
sub get_granular_array{
    my $ontology_code = shift;
    my $row = shift;
    my $evidence_exclusion_1 = shift;
    my $evidence_exclusion_2 = shift;
    my $parents_ref = shift;
    my @columns = split(/\,/, $row);
    my $size = @columns;
    my $gene_id  = $columns[0];
    my $gene_name= $columns[1];
    my @parent_list=();
        foreach my $col (@columns){
             next if ($col =~m/\[$evidence_exclusion_1\]/);
             next if ($col =~m/\[$evidence_exclusion_2\]/);
             next if ($col !~ /$ontology_code/);
             chomp($col);
             my $id = $col;
             my $evidence ="";
             ($evidence) = $id =~ /\[(.*)\]/;
             if ($evidence) {
                 $id =~ s/\[$evidence\]//g;
              }
             $id =~ s/^\s+//;
             $id =~ s/\s+$//;
             my $parent_array_ref = get_parents($id, $parents_ref);
             my @parent_array = @{$parent_array_ref};
              foreach my $parent (@parent_array){
                chomp($parent);
                $parent =~ s/^\s+//;
                $parent =~ s/\s+$//;
                push(@parent_list, $parent);
              }           
            }
            my @unique_parents = uniq(@parent_list);
             foreach my $p (@unique_parents){
              my $pointer = 0;
               while ($pointer <= $#columns) {
                   my $term = $columns[$pointer];
                   my ($evidence) = $term =~ /\[(.*)\]/;
                   if ($evidence) {
                       $term =~ s/\[$evidence\]//g;
                    }
                   $term =~ s/^\s+//;
                   $term =~ s/\s+$//;
                   if ($term eq $p) {
                       splice(@columns, $pointer, 1);
                   }
                   else {
                       $pointer++;
                   }
                 }
              }
   return \@columns;    
}
#
# get_parents returns a reference to an array of parents given an ontology identification code and 
#             given a hash containing a comma separated string of identification codes of parents 
#
sub get_parents{
my $term = shift;
my $parent_hash = shift;
my $comma = "\,";
my %parents = %$parent_hash;
my @parent_array=();
my $target_string = $parents{$term};
if ($target_string){
 my @targets = split($comma, $target_string);
 foreach my $target (@targets){
  push(@parent_array, $target);
 }
}
my @unique_parent_array = uniq(@parent_array);
return \@unique_parent_array;
}
#
# get_ontology_term returns a term given an identification code and reference to a hash of ontology terms.
#
sub get_ontology_term{
 my $ID =  shift;
 my $ontology = shift;
 my $term = $ontology->{$ID};
 return $term;
}
#
# get_ontology_parents_children returns references to hashes of ontology parents and children given the 
#                               input file (an ACE file) and ontology (e.g., "GO"). 
#
sub get_ontology_parents_children {
my $input_file = shift;
my $ontology = shift;
my %parents=();
my %children=();

if ($ontology){
 if (-e $input_file){
  my $ace_file = read_file($input_file);
  my (@terms) = split/$ontology\_term\:/, $ace_file;
  foreach my $term (@terms){
  
  my $wb="";
  my $definition="";
  $term =~ s/^\s+//;
  $term =~ s/\s+$//;
  if ($term =~ m/^\"$ontology\:(.*)/) {
      $wb = $1; 
      $wb=~s/^\s+//; 
      $wb=~s/\s+$//; 
      $wb=~s/\"//g;
      chomp($wb); 
      $definition = "$ontology\:" . $wb;
  }
           
  my (@lines) = split/\n/, $term;
  foreach my $line (@lines) {
   if ($line =~/Descendent\t(.*)/){
      my $descendent = $1; 
      $descendent=~s/^\s+//; 
      $descendent=~s/\s+$//; 
      $descendent=~s/\"//g;
      chomp($descendent);   
      if ($definition){
          $children{$definition} .= $descendent . "\,"; 
      }   
    }
   if ($line =~/Ancestor\t(.*)/){
      my $ancestor = $1;
      $ancestor=~s/^\s+//; 
      $ancestor=~s/\s+$//; 
      $ancestor=~s/\"//g;
      chomp($ancestor);  
      if ($definition){
          $parents{$definition} .= $ancestor . "\,";
      }    
    }
   }
  } # foreach term
 }
}
 return (\%parents, \%children);

}
#
# get_array requires an $ontology_code (e.g., "GO"), reference to the array that must be filled or augmented
#            an input file (parsed gene association file), optional evidence exclusions and 
#           returns a reference to an array of ontology terms.
#           This is useful to obtain a list of ontology terms used in a file. 
#
sub get_array{
    my $ontology_code = shift;
    my $array_ref = shift;
    my $input_file = shift;
    my $evidence_exclusion_1 = shift;
    my $evidence_exclusion_2 = shift;
    my @array=@{$array_ref};
    if (-e $input_file){
      my @rows = read_file($input_file);
      foreach my $row (@rows){
       my @columns = split(/\,/, $row);
       my $size = @columns;
       my $gene_id  = $columns[0];
       my $gene_name= $columns[1];
        foreach my $col (@columns){
             next if ($col =~m/\[$evidence_exclusion_1\]/);
             next if ($col =~m/\[$evidence_exclusion_2\]/);
             next if ($col !~ /$ontology_code/);
             chomp($col);
             my $id = $col;
             my $evidence ="";
             ($evidence) = $id =~ /\[(.*)\]/;
             if ($evidence) {
                 $id =~ s/\[$evidence\]//g;
              }
             push(@array, $id);
        }
      }
     }
   return \@array;    
}
#
# get_ontology_hash returns a reference to hash of ontology terms for an identification code look-up.
#                   requires a file containing the identification codes and their terms.
#                   $ontology_code is "GO" for gene ontology.
#                   One may need an optional list of alternate identification codes ($altid_file).
#
sub get_ontology_hash{
my $ontology_code = shift;
my $file  = shift;
my $altid_file = shift;
my %ontology;
chomp($ontology_code);
$ontology_code =~ s/^\s+//;
$ontology_code =~ s/\s+$//;
if (-e $file){
  my @lines = read_file($file); 
  foreach my $line (@lines){
  chomp($line);
  my $key   = $line;
  my $value = $line;
  $key  =~ /\(($ontology_code.+?)\)/;
  $key = $1; 
  $value =~ s/\($ontology_code.*//g;
  chomp($value);
  $value =~ s/^\s+//;
  $value =~ s/\s+$//;
  $ontology{$key} = $value;
 }
}
#
if (-e $altid_file){
   my @altid_lines = read_file($altid_file);
   foreach my $row (@altid_lines){
       chomp($row);
       my ($goid, $alt, $name) = split(/\t/, $row);
       my @alts = split(/\-/,$alt);
       foreach my $alt_id (@alts){
        chomp($name);
        chomp($alt_id);
        $alt_id =~ s/^\s+//;
        $alt_id =~ s/\s+$//;
        $name =~ s/^\s+//;
        $name =~ s/\s+$//;
        $ontology{$alt_id} = $name;
      }
 }
}
return \%ontology;
}
sub get_species_ensp_hash{
   my $species = shift;
   my $file  = shift;

  my $csv = Text::CSV->new ({
     binary    => 1,
     auto_diag => 1,
     sep_char  => ','    # not really needed as this is the default
  });

   my %hash=();

   open(my $data, '<:encoding(utf8)', $file) or die "Could not open '$file' $!\n";
   while (my $fields = $csv->getline( $data )) {

        my $protein_id = $fields->[0]; 
        my $ensp_id = $fields->[3];
        $protein_id =~ s/^\s+//;
        $protein_id =~ s/\s+$//;
        $ensp_id =~ s/^\s+//;
        $ensp_id =~ s/\s+$//;
        $hash{$protein_id} = $ensp_id;
    }
  close($data); 
  return \%hash; 
}
sub get_gene_protein_hash{
 my $species = shift;
 my $xrefs   = shift;
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
