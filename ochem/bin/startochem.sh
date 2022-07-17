#!/bin/sh

cd /etc/ochem ; unzip -o /home/ochem.zip
cp /ochem/cfg/* .

unzip -o /ochem/ochem.release -d /ochem/tmp/cs_release/
#mv /ochem/tmp/cs_release/ochem.war /etc/ochem/ochem-tomcat/webapps/ROOT.war
#mv /ochem/tmp/cs_release/metaserver.war /etc/ochem/metaserver-tomcat/webapps/metaserver.war

echo export OCHEMEMORY=4096
echo export METAMEMORY=2048

#OCHEM
mkdir -p /etc/ochem/ochem-tomcat/webapps/ROOT
cd /etc/ochem/ochem-tomcat/webapps/ROOT
jar -xvf /ochem/tmp/cs_release/ochem.war
cd js
rm -rf lib
ln -s /etc/source/javascripts lib

#METASERVER
mkdir -p /etc/ochem/metaserver-tomcat/webapps/metaserver
cd /etc/ochem/metaserver-tomcat/webapps/metaserver
jar -xvf /ochem/tmp/cs_release/metaserver.war
cd js
ln -sf /etc/source/javascripts/jquery.flot.js .
ln -sf /etc/source/javascripts/jquery-1.4.1.min.js .

sh /etc/ochem/ochem/bin/metaserver-tomcat start
sh /etc/ochem/ochem/bin/ochem-tomcat start
