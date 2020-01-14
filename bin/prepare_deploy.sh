#!/bin/bash
set -x
# remove all the python cruft
find . | grep -E "(__pycache__|\.pyc|\.pyo$)" | xargs rm -rf
mkdir deploy
cp -R notebooks/* deploy/
ls -R deploy
