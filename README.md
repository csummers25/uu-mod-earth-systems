# Modelling Earth Systems - Project Code
A python reworking of the visco-elasto-plastic geodynamics code developed in the book *Numerical Geodynamic Modelling* by Taras Gerya [1]. It uses a finite difference staggered grid approach to solve the Stoke's, continuity and temperature equations, along with markers for advection.

## Installation
If you haven't already, made sure you have python 3 installed (you can do this by installing Anaconda from https://www.anaconda.com/download, which will also include the packages required for this code).

The required python packages are listed in the `environment.yml` file.  To create a new conda environment with these packages installed, open a terminal (or on Windows the `anaconda_prompt` app from within the Anaconda naviagator), navigate to the directory in which you have downloaded this repository using `cd ./path/to/directory/on/your/computer/uu-mod-earth-systems` and then run the command `conda env create -f environment.yml` in the top directory of the repo.  This will create an environment called `modearth`, which can then be activated using `conda activate modearth` in the terminal (or by selecting it in the environments tab in Anaconda Navigator).  

## Running different model setups
The file `main.py` contains the central code, to run a model you can simply run this file from your choice of python IDE/interpreter.  The environment installed in the previous step includes Spyder so we recommend using this, after activating your environment as described in the install instructions you can launch spyder by typing `spyder` in the terminal (or by lauching it through the Anaconda Navigator).  You can then open the `main.py` file in spyder and run it using the run icon in the top bar. 

Various model setups can be found in the `models` directory.  In each sub-directory you will find three files: `setup.py`, `material_properties.txt` and `visualisation.py`.  These contain the code and required material parameters for that model, in a function called `initializeModel`, which is called by `main.py` to setup the simulation. The first two files sets the initial conditions, boundary conditions and parameter values for the simulation. The last is used for the visualisation of the produced output. 

To run a given setup, change the "model_name" at the top in `main.py` to point to the directory of the model you wish to run. For example, to run the model in the directory `Subduction`, we change the model_name to "Subduction"

## Changing model parameters
In the exercises you will be asked to change parameters for the different model setups. The best way to do this and keep everything organised is: Choose the sub-directory of the setup you want to use. Copy the entire folder and rename it to a clear, descriptive name. In the new sub-directory you can modify `setup.py` and `material_properties.txt` for your own implementation. In the main script, change the variable "model_name" to match the name of your new sub-directory. Now you can run the model, and it will use your modified setup. 

## References
[1] Gerya T. MATLAB program examples. In: Introduction to Numerical Geodynamic Modelling. Cambridge University Press; 2019:425-437. DOI: https://doi.org/10.1017/9781316534243.024 
