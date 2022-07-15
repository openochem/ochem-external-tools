# 1) Copy the udf module shared library mysql_udfmodule.so to the Mysql plugin directory
#   (usually something like /usr/lib64/mysql/plugin)
# 2) Verify all additional shared libraries can be found (run ldd on the module shared library, as mysql user)
#    If your SSL system libs are outdated, you may heve to copy the libssl.so, libcrpyto.so and libsasl2.so
#    libs from the toolkit package into the plugin directory to satisfy all dependencies
# 3) make sure AppArmor is not interfering (see http://forums.mysql.com/read.php?117,298997,300655#msg-300655)
# 4) Finally, run below SQL statements to define the additional module functions.
#    Run as "mysql <mysqludfinit.sql" or use "source mysqludfinit.sql" on the mysql command line
#    You probably need DB admin permissions for that. 

use mysql; 

#delete from mysql.func;

# Set up function definitions
drop function if exists ens_string_property;
drop function if exists ens_blob_property;
drop function if exists ens_double_property;
drop function if exists ens_long_property;

drop function if exists reaction_string_property;
drop function if exists reaction_blob_property;
drop function if exists reaction_double_property;
drop function if exists reaction_long_property;

drop function if exists match_formula;
drop function if exists match_substructure;
drop function if exists match_substructure_ex;
drop function if exists match_superstructure;
drop function if exists similarity;

drop function if exists string_script;
drop function if exists double_script;
drop function if exists long_script;
drop function if exists void_script;
drop function if exists blob_script;

drop function if exists string_aggregate_script;
drop function if exists blob_aggregate_script;
drop function if exists double_aggregate_script;
drop function if exists long_aggregate_script;
drop function if exists void_aggregate_script;

drop function if exists print_bitvector;

create function ens_string_property returns string soname 'mysql_udfcactvsmodule.so';
create function ens_blob_property returns string soname 'mysql_udfcactvsmodule.so';
create function ens_double_property returns real soname 'mysql_udfcactvsmodule.so';
create function ens_long_property returns integer soname 'mysql_udfcactvsmodule.so';

create function reaction_string_property returns string soname 'mysql_udfcactvsmodule.so';
create function reaction_blob_property returns string soname 'mysql_udfcactvsmodule.so';
create function reaction_double_property returns real soname 'mysql_udfcactvsmodule.so';
create function reaction_long_property returns integer soname 'mysql_udfcactvsmodule.so';

create function match_formula returns integer soname 'mysql_udfcactvsmodule.so';
create function match_substructure returns integer soname 'mysql_udfcactvsmodule.so';
create function match_substructure_ex returns integer soname 'mysql_udfcactvsmodule.so';
create function match_superstructure returns integer soname 'mysql_udfcactvsmodule.so';
create function similarity returns integer soname 'mysql_udfcactvsmodule.so';

create function string_script returns string soname 'mysql_udfcactvsmodule.so';
create function blob_script returns string soname 'mysql_udfcactvsmodule.so';
create function long_script returns integer soname 'mysql_udfcactvsmodule.so';
create function double_script returns real soname 'mysql_udfcactvsmodule.so';
create function void_script returns integer soname 'mysql_udfcactvsmodule.so';

create function print_bitvector returns string soname 'mysql_udfcactvsmodule.so';

create aggregate function string_aggregate_script returns string soname 'mysql_udfcactvsmodule.so';
create aggregate function blob_aggregate_script returns string soname 'mysql_udfcactvsmodule.so';
create aggregate function long_aggregate_script returns integer soname 'mysql_udfcactvsmodule.so';
create aggregate function double_aggregate_script returns real soname 'mysql_udfcactvsmodule.so';
create aggregate function void_aggregate_script returns integer soname 'mysql_udfcactvsmodule.so';

# Show the currently defined set of functions
# select * from mysql.func;

