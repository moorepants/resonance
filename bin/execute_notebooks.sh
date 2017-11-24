#!/bin/bash -xe

# %matplotlib notebook will not cause plots to render in html output, so swap it.
if [ "$TRAVIS" = "true" ]
then
	echo "hello"
	sed -i -- 's/%matplotlib notebook/%matplotlib inline/g' notebooks/*.ipynb
fi

for f in notebooks/*.ipynb; do
  jupyter nbconvert --debug --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 --to=html "$f"
done
