### Description ###

Set of scripts for performing several types of MD calculations for polymer physics using [openMM](http://openmm.org/) software.

### Setup ###

Set virtual env for python 3.7.x using [pyenv](https://github.com/pyenv/pyenv-installer), with separate virtual env.
Compile openMM from sources: [How-to](http://docs.openmm.org/latest/userguide/library.html#compiling-openmm-from-source-code)

Probably you will need to correct path to python executable, and instalation localization at ccmake.

Clone repo:

    git clone git@bitbucket.org:4dnucleome/md_soft.git

Install dependencies:

    pip install -r requirements.txt

### Running simulations
Assuming your config file is named `config.ini`

    run.py -c config.ini

### Other topics

##### Preparing initial structure
Your simulation need initial structure. You may find useful scripts for generating initial structures in [this repository](https://bitbucket.org/mkadlof/structuregenerator).

##### General tips
Distance Constraints on consecutive beads doesnt work with harmonic flat angle. Use harmonic bond instead. 

#### Contact
Michał Kadlof <m.kadlof@cent.uw.edu.pl>