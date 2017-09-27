#!/bin/bash
set -x  # print all executed commands

# copy ipynb source files to notebooks dir
for srcpath in notebook-src/*.ipynb
do
	filename=$(basename "$srcpath")
    cp "$srcpath" "notebooks/$filename"
done

# copy figures over
cp notebook-src/fig/* notebooks/fig/

# build rst files into ipynb
for srcpath in notebook-src/*.rst
do
	filename=$(basename "$srcpath")
	rst2ipynb "$srcpath" "notebooks/${filename%.*}.ipynb"
done

# %matplotlib notebook will not cause plots to render, so swap it.
if [ "$TRAVIS" = "true" ]
then
	echo "hello"
	sed -i -- 's/%matplotlib notebook/%matplotlib inline/g' notebooks/*.ipynb
fi
