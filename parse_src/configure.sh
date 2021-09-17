#!/bin/bash


pip install flask requests numpy Cython pyOpenSSL pandas
#pip install cyvcf2
cd cmd2web/

git submodule init
git submodule sync
git submodule update --recursive

cd ..

rm -rf cmd2web/src/web_src
rm -rf cmd2web/src/Parse.f


cp -r web_src cmd2web/src/
cp -r parse_src/Parse.f cmd2web/src/web_src/static/js/
