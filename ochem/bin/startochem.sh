#!/bin/sh

cd /etc/ochem ; unzip -o /home/ochem.zip
cp /ochem/cfg/* .

unzip -o /ochem/ochem.release -d /ochem/tmp/cs_release/
mv /ochem/tmp/cs_release/ochem.war /etc/ochem/ochem-tomcat/webapps/ROOT.war
mv /ochem/tmp/cs_release/metaserver.war /etc/ochem/metaserver-tomcat/webapps/metaserver.war

echo export OCHEMEMORY=4096
echo export METAMEMORY=2048

sh /etc/ochem/ochem/bin/metaserver-tomcat start
sh /etc/ochem/ochem/bin/ochem-tomcat start
