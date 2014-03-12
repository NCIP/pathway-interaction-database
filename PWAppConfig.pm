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

  $ENV{'ORACLE_HOME'} = "/oracle/product";

## Global variables for database schema names

use constant DB_USER         => "useid";
use constant DB_PASS         => "password";
use constant DB_INST         => "db";
use constant SCHEMA          => "pid";

use constant DOT_PGM         => "/usr/local/bin/dot";
use constant PW_SERVER_PORT         => "8001";

use constant INIT_DATA_HOME          => "/PID/data/";

use constant CACHE_ROOT              => "/PID/data/cache/";

use constant PW_CACHE_PREFIX        => "PW";

use constant PW_COMMENT_LOG         => "/PID/data/cache/pw_comment.log";

# Perl equivalents of Zope dtml-var's BASE and IMG_DIR
use constant BASE     => "";
use constant IMG_DIR  => "/images";

use constant SAFE_IPS =>
 "127.0.0.1," .           ## localhost
   "xxx.xxx.xxx.xxx" ;        ## 
  

######################################################################
1;
