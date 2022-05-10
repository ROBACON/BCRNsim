
# BCRN simulations using COPASI in Python

This package allows for streamlined simulations of Biochemical reaction networks (BCRNs) using Copasi (https://copasi.org)
solvers.
Further, merging several stochastic simulations, parallelization of runs, and plotting are supported. 

Organization:
1. Requirements
2. Directories
3. Basic Usage
4. Example Models


## 1. Requirements

Copasi (https://copasi.org) must be installed and added to the PATH environment variable, such that CopasiSE can be called via the command line from anywhere.

Python3 and the relevant python libraries must be installed. One can install all or most of the required python libraries via "pip3 install -r requirements.txt." Any remaining libraries must be installed manually using pip following the relevant error message.


## 2. Directories

- `simulation`: the simulation framework 
  - `params`: parameter files
  - `models`: set of available models and a builder to translate them into Copasi
  - `lib`: utility scripts, suchs as statistics and parsing
    - `cps` and `output`: internal temporary directories to store Copasi cps files and output
  - `figure_styles`: styling of figures
- `data`: converted data file from the work "A synthetic communication system uncovers self-jamming of bacteriophage transmission" by Amit Pathania, Corbin Hopper, Amir Pandi, Matthias Függer, Thomas Nowak, Manish Kushwaha.


## 3. Basic Usage

Running a single simulation: in the `simulation` run
	`python3 run.py params/PARAM_FILE`

With the default parameter file `params/twophase.txt`, the simulation should run and the image `Two-phase model.png` should be generated in the current (`simulation`) directory. 

Modify the chosen PARAM_FILE to rapidly pick a model and specify simulation parameters. For instance, one can change the algorithm used for simulation, initial concentrations of species, or kinetic rates of the model. 

In the case of the `params/twophase.txt` example, the model parameters are found in `params/twophase.txt`, the model, i.e, the BCRN reactions, species, and events, are generated by `models/twophase.py`, and the plotting and output parameters are specified in `params/twophase-plot.json`.

For a more advanced example of a model that maps several sub-species to abstract species, as well as more refined plotting and output options, see the example `params/Recipient_Kan.txt`.
This example also includes external data from experiments in the output plot.

## 4. Example Models

The following models have been used in the work "A synthetic communication system uncovers self-jamming of bacteriophage transmission" by Amit Pathania, Corbin Hopper, Amir Pandi, Matthias Függer, Thomas Nowak, Manish Kushwaha.


### 4.1. E. coli growth models

Zero-resource, one-resource, and two-resource models
for growth of E. coli (Supp. Fig. Box-1 in the paper)
can be simulated via

`python3 run.py params/zerophase.txt` 
`python3 run.py params/onephase.txt`
`python3 run.py params/twophase.txt`

Output plots are generated in the current directory.


### 4.2. Donor E. coli model

The phage secretion model (Fig. 1e in the paper) is simulated via

`python3 run.py params/Donor_cfu.txt`


### 4.3. Receiver E. coli model

The model for Receiver cells being infected by phages (Fig. 2 in the paper) is simulated via

`python3 run.py params/Recipient_Kan.txt`


### 4.4. Combined Donor and Receiver model

The combined model that contains Donors secreting phages
that infect Receivers (Fig. 3 in the paper) is simulated via

`python3 run.py params/DonorRecipient.txt`


### 4.5. Long-term infection with M13 phages

The model of M13 infecting E. coli in a gut-like environment (Fig. 5 in the paper)
is simulated via

`python3 run.py params/Wildtype.txt`


