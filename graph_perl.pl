#!/usr/bin/env perl
use warnings;
use strict;
use Devel::Graph;
use File::Slurp;
use Graph::Easy;

my @graphics_array = qw(bmp canon dot gv xdot xdot1.2 xdot1.4 cgimage 	
                      cmap eps exr fig gd gd2 gif gtk ico imap 
                      cmapx imap_np cmapx_np ismap jp2 jpg jpeg jpe 	
                      pct pict pdf pic plain plain-ext png pov ps 
                      ps2 psd sgi svg svgz tga tif tiff tk vml 
                      vmlz vrml wbmp webp xlib x11);

    my $argc = @ARGV;
    my $TXT= ".txt";
    my $HTML = ".html";

    if ($argc < 3) {
     print "\nUsage: graph_perl.pl file_name output_format output_directory\n";
     exit;
   }
        
        my $file = $ARGV[0];
        my $graphics_ext = lc $ARGV[1];
        my $dir = $ARGV[2] . "\/" ;
     
     my $found = 0;
     foreach my $format (@graphics_array){   
       if ($format =~/$graphics_ext/i){
         $found = 1;
       }
     }
     if ($found == 0){
      print "\nOutput format not found\.  Use one of the following output formats\n";
      foreach my $format (@graphics_array){
        $format =~s/\"//g;
        $format =~s/\,//g;
        print "$format  "; 
       }
      print "\n\nUsage: graph_perl.pl file_name output_format output_directory\n\n";
     exit;
     } 
        my $name = $file;
        $name =~ s{.*/}{};      # removes path  
        $name =~ s{\.[^.]+$}{}; # removes extension
        my ($ext) = $file =~/(\.[^.]+)$/;
        if ( ($ext=~/pl/i) or ($ext=~/pm/i) ){

        my $grapher = Devel::Graph->new();
        my $output_file = $file;
        my $output_text_file = $file;
        my $output_html_file = $file;

        my $replace = $ext;
           $output_file = $dir . $name . "\." . $graphics_ext;
           $output_text_file = $dir . $name . $TXT;
           $output_html_file = $dir . $name . $HTML;
        print "file to graph is $file\.\n\n";
        print "output files are in directory\, $dir\:\n $output_text_file\n $output_html_file\n $output_file\n\n"; 
        my $graph = $grapher->decompose($file);
           $graph->catch_warnings(1);
        my @warnings = $graph->warnings();
           if (@warnings) {
             print "There were some warnings\.\n\n";
           }
        my $out_text = $graph->as_ascii();
           $out_text .= "\n";
        write_file($output_text_file, $out_text);

        my $graph2 = $grapher->as_graph();
           $graph2->catch_warnings(1);
        my @warnings2 = $graph2->warnings();
           if (@warnings2) {
             print "There were some warnings\.\n\n";
           }
         $graph2->set_attribute('node.if', 'fill', 'yellow');    # if blocks:
         $graph2->set_attribute('node.for', 'fill', 'cyan');  # for blocks: 
         $graph2->set_attribute('edge.true', 'style', 'bold');# true edges: 
        my $flow = Devel::Graph->new();
           $flow->decompose( $file );
           $flow->catch_warnings(1);
        my @warnings3 = $flow->warnings();
           if (@warnings3) {
             print "There were some warnings\.\n\n";
           }

        my $chart = $flow->as_flowchart();
        my $html = $chart->as_html_file();
         $html=~s/\<html\>/\<html\>\<title\>$name\<\/title\>/g;
         write_file($output_html_file, {binmode => ':utf8'}, $html);
        
        my $graphviz = $graph2->as_graphviz();
        my $DOT;
        open $DOT, "|dot -T$graphics_ext -o $output_file" or die ("Cannot open pipe to dot: $!");
        print $DOT $graphviz;
        close $DOT;
        
        $flow->finish();

        } else {
          print "use is only for perl scripts at this time";
        }
exit 0;
