#!/bin/bash

for srcpath in notebook-src/*.rst
do
	filename=$(basename "$srcpath")
	rst2ipynb $srcpath "notebooks/${filename%.*}.ipynb"
done

cp notebook-src/fig/* notebooks/fig/
