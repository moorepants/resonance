#!/bin/bash -xe

# %matplotlib notebook will not cause plots to render in html output, so swap it.
if [ "$TRAVIS" = "true" ]
then
	echo "hello"
	find notebooks/ -name "*.ipynb" -exec sed -i -- 's/%matplotlib notebook/%matplotlib inline/g' {} \;
fi

for dir in notebooks/ notebooks/*/ ; do
	if [ $dir != "notebooks/scratch/" ]; then
		echo "Converting in directory $dir."
		echo "============================="
		cd $dir
		for nb in *.ipynb ; do
			echo "Converting $nb."
			echo "==============="
			jupyter nbconvert --debug --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 --to=html "$nb"
		done
		cd -
	fi
done
