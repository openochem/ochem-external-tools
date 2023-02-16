#!/bin/sh

for var in "$@"
do
    if [[ $var =~ ^ochemenv(.+) ]]; then
        port=`echo "$var" | sed 's/ochemenv//'`
        port=${port/\/}
    else
        export "$var"
    fi
done

cd /etc/ochem ; unzip -o /home/ochem.zip
cp /ochem/cfg/* .

unzip -o /ochem/ochem.release -d /ochem/tmp/cs_release/

if [ -z "$port" ]; then 
		echo "using default port 8080"; 
	else 
		echo "using port '$port'";
		sed -i -r "s|8080|$port|g" /etc/ochem/ochem-tomcat/conf/server.xml
		summa=`expr 9009 + $port - 8080`
		echo "using port '$summa'";
		sed -i -r "s|8009|$summa|g" /etc/ochem/ochem-tomcat/conf/server.xml
		summa=`expr 9006 + $port - 8080`
		sed -i -r "s|8005|$summa|g" /etc/ochem/ochem-tomcat/conf/server.xml
fi

#OCHEM
mkdir -p /etc/ochem/ochem-tomcat/webapps/ROOT
cd /etc/ochem/ochem-tomcat/webapps/ROOT
jar -xvf /ochem/tmp/cs_release/ochem.war >/dev/null 2>/dev/null

cd js
rm -rf lib; mkdir lib
cp -r /etc/source/javascripts/* lib
rm -rf lib/jquery-1.4.1.min.js
cp /etc/source/cdk/* /etc/ochem/ochem-tomcat/webapps/ROOT/WEB-INF/lib/ # can be used to update CDK
cp /etc/source/chem/* /etc/ochem/ochem-tomcat/webapps/ROOT/WEB-INF/lib/ # any additional code to be copied can be added here

if [ -z "$port" ]; then 
	#METASERVER
	mkdir -p /etc/ochem/metaserver-tomcat/webapps/metaserver
	cd /etc/ochem/metaserver-tomcat/webapps/metaserver
	jar -xvf /ochem/tmp/cs_release/metaserver.war >/dev/null 2>/dev/null
	cd js
	cp /etc/source/javascripts/jquery.flot.js .
	cp /etc/source/javascripts/jquery-1.4.1.min.js .
	sh /etc/ochem/ochem/bin/metaserver-tomcat start
fi

sh /etc/ochem/ochem/bin/ochem-tomcat start
