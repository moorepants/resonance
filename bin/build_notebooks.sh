#!/bin/bash
set -x  # print all executed commands
for srcpath in notebook-src/*.rst
do
	filename=$(basename "$srcpath")
	rst2ipynb "$srcpath" "notebooks/${filename%.*}.ipynb"
done
# %matplotlib notebook will not cause plots to render, so swap it.
if [ "$TRAVIS" = "true" ]
then
	sed -i -- 's/%matplotlib notebook/%matplotlib inline/g' "notebooks/*.ipynb"
fi
cp notebook-src/fig/* notebooks/fig/
