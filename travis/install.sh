#!/bin/bash

if [ $TRAVIS_OS_NAME = 'osx' ]; then
    sudo pip install --upgrade pip
    pip install --user numpy
    pip install --user pygmo
else
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh;
    bash miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    hash -r
    conda config --set always_yes yes --set changeps1 no
    conda update -q conda
    conda info -a
    conda create -q -n test-env python=$TRAVIS_PYTHON_VERSION
    source activate test-env
    conda install -c anaconda numpy
    conda install -c conda-forge pygmo
fi
