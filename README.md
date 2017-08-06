# wear_modelling

Geomechanical boundary element modelling of the behavor of 2 dimensional profiles in subsect to stresses. The model is built to assess the behavior of faults and, specifically, the failure  of geometrical asperities on slip surfaces. Matlab sripts develloped to work with Fric2d and GROW projects develloped by Michele Cooke and Jess McBeth. Fric2d computes tractions and displacement along the fault profile. GROW computes the propagations of cracks by incrementally adding fault elements so as to minimize work. The matlab script was built to provive simple input structures, runs scripts and extract output information, assess the tendency and spatial locale for failure and graphically display all the information.

## Getting Started

Dowload a copy of the master branch and unzip it in desired directory. Open matlab (2016a or later).

### Prerequisites

UNIX or MAC system. Does not work with windows systems. Matlab 2016a or later (https://www.mathworks.com/). Perl (should be installed on all UNIX or MAC systems)

### Installing

Once matlab is open, you should be able to run the command:

```
sample_run
```
This will automatically compile a working version of the fric2d source code if it is not in the working directory

Depending of the model set up, an output figure should appear after ~ a minute or so

If a model ran to completion the output will be stored in a directory named 'output', important files will be the input (the shortest file with extension .in) file which is the input file to fric2d and 'grow_specs.in' which is the input to GROW. The figure will be stored in a .fig file, the will also be a matlab workspace saved with input paramaters that describe the material properties and boundary conditions used for the run.

### Changing input parameters
The work flow is essentially composed of two matlab files. 'fric2d_workflow' takes as input x coordinates, y coordinates and a fileName. These are supplied by 'sample_run.m'. In this file you may change the input geometry. It also has the ability of sending out batch commands to fric2d_workflow. To change the number of runs change the line:

```
numberOfRuns = 1;
```
if doing multiple runs you may iterate the height, length, and shape of asperities by uncommenting one or many of:

```
% absoluteAsperityLenghtArray  = 0.01./(1:numberOfRuns);
% absoluteAsperityHeightArray = 0.005./(1:numberOfRuns);
% asperityTypeArray            = {'sine','triangle', 'box','step','jog'}; if length(asperityType ~= numberOfRuns); error(asperityTypes must be of length number of runs
```
To change material properties or more specific behavioral parameters of fric2d or GROW open the 'fric2d_workflow.m' file in matlab. To simplify the interface press:

```
ctl +=
```
to collapse all the function. Expand the function

```
function [inputParameters, GROWInputParameters] = userInput()
    inputParameters = [];
```
Here is should be evident how to change material propeties and boundary conditions. 

Keep in mind that there is a complex interplay between boundary conditions, elemeent size, material properties and displacement on endbits.

## Authors

**Kelian Dascher-Cousineau** workflow, **Michele Cooke** - Fric2d, **Jess Mecbeth** GROW

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
