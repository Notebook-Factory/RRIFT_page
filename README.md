
# RRIFT_notebooks

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Notebook-Factory/RRIFT_notebooks/HEAD?filepath=RRIFT.ipynb)

Paper: https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.27913 <br> 
Code: https://github.com/MPUmri/RRIFT <br> 
Note: MATLAB. Reproduces ALL figures

# Setup instructions

## Prerequisites

Make sure before running the notebook you download ```simResults.mat```, ```simResultsTRes15-varKtRR.mat```, ```simResultsTRes15-varVeRR.mat``` and ```simMap.mat``` from the [OSF repository](https://osf.io/wr3kf/files/) and place them inside ```./RRIFT/data/```. Also make sure you have [```Anaconda 3```](https://www.anaconda.com/products/individual) installed on your system as well as ```MATLAB R2020b```.

## Installing MATLAB

You can download MATLAB from [here](https://uk.mathworks.com/downloads/).

## Making a conda environment

* ```conda create -n rrift_notebook python=3.6``` (you can set ```rrift_notebook``` to whatever you want the name of your environment to be)
* ```conda activate rrift_notebook```

## Installing the MATLAB kernel for Jupyter

* ```pip install matlab_kernel```

To check if the kernel is properly installed use the command ```jupyter kernelspec list``` which should list both MATLAB and Python if installed correctly.

## MATLAB-side configuration 

The MATLAB executable needs to be exposed to Jypiter.
You need to go to the directory where MATLAB is installed and navigate inside `extern/engines/python`. 
For example, the path should look something like this on Windows: `M:\MATLAB\extern\engines\python` or on Linux: `/usr/local/MATLAB/R2020b/extern/engines/python`. 

Once you've navigated to that location you need to install the Python engine with the command:

* ```python setup.py install```

## Setting up SoS and MATLAB 

* ```conda install sos sos-pbs -c conda-forge```
* ```conda install sos-notebook jupyterlab-sos sos-papermill -c conda-forge```
* ```conda install sos-python sos-matlab -c conda-forge```

## Installing Python packages

* ```pip install plotly=4.14.3 numpy=1.19.5 scipy=1.6.0 statsmodels=0.12.1 pandas=1.2.0```

## Start a Jupyter session

* ```jupyter notebook```
