#!/bin/bash -xe

# %matplotlib notebook and widget will not cause plots to render in html output, so swap it.
if [ "$TRAVIS" = "true" ]
then
	find notebooks/ -name "*.ipynb" -exec sed -i -- 's/%matplotlib notebook/%matplotlib inline/g' {} \;
	find notebooks/ -name "*.ipynb" -exec sed -i -- 's/%matplotlib widget/%matplotlib inline/g' {} \;
	sed -i 's/.ipynb/.html/g' notebooks/index.ipynb
fi

cd notebooks
for dir in */ ; do
	zip -r ${dir:0:-1}.zip $dir
done
cd -

for dir in notebooks/*/ notebooks/ ; do
	if [ $dir != "notebooks/scratch/" ]; then
		echo "Converting in directory $dir."
		echo "============================="
		cd $dir
		for nb in *.ipynb ; do
			echo "Converting $nb."
			echo "==============="
			jupyter nbconvert --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 --to=html "$nb"
		done
		if [ -d ".ipynb_checkpoints" ] ; then
			rm -r .ipynb_checkpoints
		fi
		cd -
	fi
done
