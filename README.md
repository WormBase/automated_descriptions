concise_descriptions
====================

automated concise descriptions

The requirements for automated concise descriptions are outlined online:
	http://wiki.wormbase.org/index.php/Generation_of_automated_descriptions

The standard operation procedure for the software is itemized online:
	http://wiki.wormbase.org/index.php/Documentation_for_workflow_and_scripts

The software is available online: 
	https://github.com/WormBase/automated_descriptions

Output can be found here: http://textpresso-dev.caltech.edu/concise_descriptions/

The perl module dependencies are:
 AcePerl, Carp, DBI, File::Slurp, LWP::Simple, LWP::UserAgent, List::MoreUtils, List::Util, OBO::Parser::OBOParser, POSIX, Spreadsheet::XLSX, Switch and Text::CSV 

Perl 5.8.x (or higher) and R 3.x are used on Ubuntu and Red Hat Enterprise Linux platforms. 

GNU Parallel 20130922 (or higher), https://www.gnu.org/software/parallel/, is required as well.

Enjoy !

R. Kishore and J. Done, California Institute of Technology, 2016
