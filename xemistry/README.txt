Required to add xemistry functionality.

mysql_udfcactvsmodule.so (version for a correct linux system can be downloaded from xemistry for non-commercial use)

https://www.xemistry.com/academic/
 
should be installed  to plugin directory of mysql/mariadb, e.g. (depening on the target architecture)

/usr/lib/x86_64-linux-gnu/libmariadb3/plugin/

It was also required to add:

/usr/lib/libidn.so.11
/usr/lib/libidn.so.11.6.16

after installation init mariadb/mysql with xemistry_init.sql

mysql < xemistry_init.sql

