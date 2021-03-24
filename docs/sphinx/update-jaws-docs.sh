#!/bin/bash
set -euo pipefail

# getting run in jaws/ root dir.
jaws_cwd=$(pwd)

# This is where we will clone the jaws-docs repo
if [[ -d "/tmp/build-docs" ]]; then
    rm -r /tmp/build-docs
fi
mkdir /tmp/build-docs && cd /tmp/build-docs

git clone git@gitlab.com:jfroula/jaws-docs.git

# copy new files 
cd jaws-docs
rsync -r $jaws_cwd/docs/sphinx/source/* docs/source/

git add . && git commit -m "updated the jaws-docs" && git push
echo "updates were successfully pushed to gitlab.com/jfroula/jaws-docs"

