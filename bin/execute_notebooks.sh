#!/bin/bash

jupyter nbextension enable exercise2/main
jupyter nbextension enable rubberband/main

for f in notebooks/*.ipynb; do
  jupyter nbconvert --debug --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 --to=html "$f"
  jupyter nbconvert --debug --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 --to=notebook "$f"
done
