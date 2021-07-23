#!/bin/bash
cd ..

git submodule init
git submodule sync
git submodule update --recursive

rm -rf cmd2web/src/web_src
rm -rf cmd2web/src/Parse.f


cp -r web_src/ cmd2web/src/
cp -r parse_src/Parse.f cmd2web/src/