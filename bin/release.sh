#!/bin/bash

version=$1
tmpd=/tmp/release_OAK_$$
host=abarth@modb.ulg.ac.be

mkdir $tmpd
wd=$PWD

cd $tmpd 
svn export svn+ssh://$host/home/svn/repos/OAK/tag/r$version OAK-$version
cd $tmpd/OAK-$version
mv config.mk.template config.mk
cd $tmpd/OAK-$version/doc
make pdf
cd $tmpd 
tar -czvf $wd/OAK-$version.tar.gz OAK-$version

scp $wd/OAK-$version.tar.gz $host:/var/lib/mediawiki/upload/Alex/OAK/release

cd $tmpd/OAK-$version/doc
mkdir html
rsync -avz . $host:/var/lib/mediawiki/upload/Alex/OAK/doc
