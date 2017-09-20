#!/bin/bash
set -x
mkdir deploy
cp -R notebooks/* deploy/
ls -R deploy
# rename `<notebook>.nbconvert.ipynb` to `<notebook>.ipynb`
for f in deploy/*.nbconvert.ipynb; do
    mv "$f" "${f%.nbconvert.ipynb}.ipynb"
done
