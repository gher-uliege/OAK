#!/bin/bash

if [ $# -ne 1 ]; then
   echo Usage $0 version_number
   exit 1
fi


version=$1
tmpd=/tmp/release_OAK_$$
host=abarth@modb.oce.ulg.ac.be

mkdir $tmpd
wd=$PWD

cd $tmpd 
svn export svn+ssh://$host/home/svn/repos/OAK/tag/r$version OAK-$version
cd $tmpd/OAK-$version
rm -Rf Python
# symbolic links
cp --dereference $(find $HOME/Assim/OAK/matlab -type l) .

mv config.mk.template config.mk
cd $tmpd/OAK-$version/doc
make pdf
cd $tmpd 
tar -czvf $wd/OAK-$version.tar.gz OAK-$version

scp $wd/OAK-$version.tar.gz $host:/var/lib/mediawiki/upload/Alex/OAK/release

cd $tmpd/OAK-$version/doc
make html
rsync -avz *pdf *html *css $host:/var/lib/mediawiki/upload/Alex/OAK/doc
