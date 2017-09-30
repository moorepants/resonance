#!/bin/bash
set -x  # print all executed commands

# copy ipynb source files to notebooks dir
for srcpath in notebook-src/*.ipynb
do
	filename=$(basename "$srcpath")
    cp "$srcpath" "notebooks/$filename"
done

# copy figures over
mkdir -p notebooks/fig/
cp -r notebook-src/fig/* notebooks/fig/

# copy data over
mkdir -p notebooks/data/
cp -r notebook-src/data/* notebooks/data/

# build rst files into ipynb
for srcpath in notebook-src/*.rst
do
	filename=$(basename "$srcpath")
	rst2ipynb "$srcpath" "notebooks/${filename%.*}.ipynb"
done

# %matplotlib notebook will not cause plots to render in html output, so swap it.
if [ "$TRAVIS" = "true" ]
then
	echo "hello"
	sed -i -- 's/%matplotlib notebook/%matplotlib inline/g' notebooks/*.ipynb
fi
