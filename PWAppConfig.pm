##
## Parameters for back-end data servers for CGAP
##

package PWAppConfig;
require Exporter;
 
@ISA     = qw(Exporter);
@EXPORT  = qw(

  DB_USER
  DB_PASS
  DB_INST
  SCHEMA

  PW_SERVER_PORT

  DOT_PGM

  INIT_DATA_HOME
  PW_CACHE_PREFIX
  CACHE_ROOT

  PW_COMMENT_LOG

  BASE
  IMG_DIR

  SAFE_IPS
);


if (-d "/app/oracle/product/dbhome/current") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/dbhome/current";
} elsif (-d "/app/oracle/product/8.1.7") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/8.1.7";
} elsif (-d "/app/oracle/product/8.1.6") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/8.1.6";
} elsif (-d "/app/oracle/product/10gClient") {
  $ENV{'ORACLE_HOME'} = "/app/oracle/product/10gClient";
}

$ENV{'PATH'} = "$ENV{'PATH'}:/app/oracle/product/10gClient:/app/oracle/product/10gClient/lib:/app/oracle/product/10gClient/instantclient:/app/oracle/product/10gClient/bin";
$ENV{'LD_LIBRARY_PATH'} = "/app/oracle/product/10gClient/lib";

## Global variables for database schema names

use constant DB_USER         => "web";
use constant DB_PASS         => "readonly";
use constant DB_INST         => "cgprod";
use constant SCHEMA          => "pid";

use constant DOT_PGM         => "/usr/local/bin/dot";
use constant PW_SERVER_PORT         => "8001";

use constant INIT_DATA_HOME          => "/share/content/PID/data/";

use constant CACHE_ROOT              => "/share/content/PID/data/cache/";

use constant PW_CACHE_PREFIX        => "PW";

use constant PW_COMMENT_LOG         => "/share/content/PID/data/cache/pw_comment.log";

# Perl equivalents of Zope dtml-var's BASE and IMG_DIR
use constant BASE     => "";
use constant IMG_DIR  => "/images";

use constant SAFE_IPS =>
  "128.231.202.171," .     ## lpgfs
  "137.187.66.221,"  .     ## lpgfs new
  "156.40.133.243,"  .     ## cbiodev24
  "137.187.66.205,"  .     ## cbiodev104
  "192.168.200.30,"  .     ## cbioapp104 from cbiodev104
  "137.187.66.207,"  .     ## cbioapp101
  "192.168.200.26,"  .     ## cbioapp101 from cbioapp101
  "137.187.66.209,"  .     ## cbioapp102
  "137.187.66.195,"  .     ## cbioapp104
  "137.187.66.222,"  .     ## ncicbstarfish
  "137.187.67.76,"   .     ## ncicbstarfish
  "137.187.67.77,"   .     ## ncicbstarfish
  "137.187.67.78," .        ## ncicbstarfish
  "137.187.67.185," .       ## cbvapp-d1016
  "137.187.182.24," .       ## cbvapp-d1017  -added if needed Bob Wysong
  "137.187.183.13" ;       ## cbiovdev5047

######################################################################
1;
