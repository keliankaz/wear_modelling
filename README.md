# wear_modelling

Geomechanical boundary element modelling of the behavor of 2 dimensional profiles in subsect to stresses. The model is built to assess the behavior of faults and, specifically, the failure  of geometrical asperities on slip surfaces. Matlab sripts develloped to work with Fric2d and GROW projects develloped by Michele Cooke and Jess McBeth. Fric2d computes tractions and displacement along the fault profile. GROW computes the propagations of cracks by incrementally adding fault elements so as to minimize work. The matlab script was built to provive simple input structures, runs scripts and extract output information, assess the tendency and spatial locale for failure and graphically display all the information.

## Getting Started

Dowload a copy of the master branch and unzip
Open matlab

### Prerequisites

Matlab 2016a or later (https://www.mathworks.com/)
Perl

### Installing

Once matlab is open, you should be able to run the command:

```
sample_run
```
This will automatically compile a working version of the fric2d source code if it is not in the working directory

Depending of the model set up, an output figure should appear after ~ a minute or so

If a model ran to completion the output will be stored in a directory named 'output', important files will be the input (the shortest file with extension .in) file which is the input file to fric2d and 'grow_specs.in' which is the input to GROW. The figure will be stored in a .fig file, the will also be a matlab workspace saved with input paramaters that describe the material properties and boundary conditions used for the run. 

## Authors

Kelian Dascher-Cousineau
Michele Cooke
Jess Mecbeth

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
