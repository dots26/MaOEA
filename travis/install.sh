#!/bin/bash

if [ $TRAVIS_OS_NAME = 'osx' ]; then
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
else
    sudo pip3 install --upgrade pip3
    pip3 install --user numpy
    pip3 install --user pygmo
fi
