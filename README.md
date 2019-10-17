### Description ###

Set of scripts for performing several types of MD calculations for polymer physics using [openMM](http://openmm.org/) software.

### Setup ###

Recommended installation is using [pyenv](https://github.com/pyenv/pyenv-installer), with separate virtual env.
This distribution was tested on `miniconda3-4.3.11`.

##### create virtual env

    cd md_soft
    pyenv install miniconda3-4.3.11
    pyenv virtualenv miniconda3-4.3.11 md_soft
    pyenv local md_soft

##### Installing dependencies

    conda install -c omnia openmm
    pip install -r requirements.txt

### Running simulations
Assuming your config file is named `config.ini`

    run.py -c config 

### Other topics

##### Preparing initial structure
Your simulation need initial structure. You may find useful scripts for generating initial structures in [this repository](https://bitbucket.org/mkadlof/structuregenerator).

#### Contact
Michał Kadlof <m.kadlof@cent.uw.edu.pl>