#!/bin/sh

for var in "$@"
do
    export "$var"
done

cd /etc/ochem ; unzip -o /home/ochem.zip
cp /ochem/cfg/* .

unzip -o /ochem/ochem.release -d /ochem/tmp/cs_release/
#mv /ochem/tmp/cs_release/ochem.war /etc/ochem/ochem-tomcat/webapps/ROOT.war
#mv /ochem/tmp/cs_release/metaserver.war /etc/ochem/metaserver-tomcat/webapps/metaserver.war

#OCHEM
mkdir -p /etc/ochem/ochem-tomcat/webapps/ROOT
cd /etc/ochem/ochem-tomcat/webapps/ROOT
jar -xvf /ochem/tmp/cs_release/ochem.war
cd js
rm -rf lib; mkdir lib
cp -r /etc/source/javascripts/* lib
rm -rf lib/jquery-1.4.1.min.js
cp /etc/source/cdk/* /etc/ochem/ochem-tomcat/webapps/ROOT/WEB-INF/lib/ # can be used to update CDK
cp /etc/source/chem/* /etc/ochem/ochem-tomcat/webapps/ROOT/WEB-INF/lib/ # any additional code to be copied can be added here


#METASERVER
mkdir -p /etc/ochem/metaserver-tomcat/webapps/metaserver
cd /etc/ochem/metaserver-tomcat/webapps/metaserver
jar -xvf /ochem/tmp/cs_release/metaserver.war
cd js
cp /etc/source/javascripts/jquery.flot.js .
cp /etc/source/javascripts/jquery-1.4.1.min.js .

sh /etc/ochem/ochem/bin/metaserver-tomcat start
sh /etc/ochem/ochem/bin/ochem-tomcat start
