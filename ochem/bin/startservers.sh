#!/bin/sh

mkdir /tmp/results
cd /tmp/
mkdir /tmp/.screen && chmod 700 /tmp/.screen
rm -rf /etc/cs_servers/*
mkdir /etc/cs_servers/server1
cd /etc/cs_servers/server1
echo "unzip /etc/ochem/release.zip"
unzip /etc/ochem/release.zip >/dev/null
cd /etc/cs_servers/
ln -s server1/tools/utils/* .
perl monitor.pl > monitor.log  &
