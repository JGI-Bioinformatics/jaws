#!/bin/bash
#MESSAGE=$1
set -euo pipefail

#if [[ !$MESSAGE ]]; then
#    echo "Usage: $0 <\"message for commit\">"
#    exit
#fi

# getting run in jaws/ root dir.
jaws_cwd=$(pwd)

if [[ -d "/tmp/build-docs" ]]; then
    rm -r /tmp/build-docs
fi
mkdir /tmp/build-docs && cd /tmp/build-docs

git clone git@gitlab.com:jfroula/jaws-docs.git

cd jaws-docs
rsync -r $jaws_cwd/docs/sphinx/source/* docs/source/

git add . && git commit -m "updated the jaws-docs" && git push
echo "updates were successfully pushed to gitlab.com/jfroula/jaws-docs"

