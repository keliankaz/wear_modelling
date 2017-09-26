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

### Under the hood:
The work flow is essentially composed of two matlab files. 'fric2d_workflow' takes as input x coordinates, y coordinates and a fileName. These are supplied by 'sample_run.m'. In this file you may change the input geometry. It also has the ability of sending out batch commands to fric2d_workflow. 

### Input parameters:

input are supplie according to matlab 'pair-wise' convention:
 ...,'numberOfRuns',1, ...      : Number of runs over which the code will be iterated
 ...,'numberOfElements',35,...  : Number of boundary elements on the fault
 ...,'asperityHeigth',0.002,... : Height of asperity in meters
 ...,'asperityLength',0.025,... : Length of asperity in meters
 ...,'shape', 'sine',...        : shape of the asperity 
                                (one of: % 'sine','triangle','box','step',
                                'jog','real profile')

Numeric inputs be iteratively changed by choosing the minimum and maximum value as the input array of length 2. I produce an array  be linear interpolation of length numberOfRuns, between these two values.

For example,

```
sample_run('numberOfRun', 5, 'asperityHeight', [0.002,0.003])
```

will run five iterations where the height is incrementally changed from the minim value, 0.002 m, to 0.003 m in heigth.

An *incorrect* usage such as

```
sample_run('numberOfRun', 5, 'asperityHeight', [0.02,0.03, 0.04])
```
wont run because the numeric specifiers must be a maximum and minimum values (i.e. of length 1 by 2).  

To change material properties or more specific behavioral parameters of fric2d or GROW, the workflow itself must be edited. Proceed with caution! Many parameters have some level of dependencies that must be taken into account when performing changes. Open the 'fric2d_workflow.m' file in matlab. To simplify the interface press:

```
ctl +=
```
to collapse all the function. Expand the function

```
function [inputParameters, GROWInputParameters] = userInput()
    inputParameters = [];
```
Here is should be evident how to change material propeties and boundary conditions. 

## Additional resources and help manuals:

Umass Geomechanical Models:
http://www.geo.umass.edu/faculty/cooke/software.html

Fric2D user manual:
http://www.geo.umass.edu/faculty/cooke/fric2d/doc.html

GROW user manual:
http://www.geo.umass.edu/faculty/cooke/grow/GROW%20User%20Manual%202015%2010-23-15.pdf

## Relevant litterature:

Cooke, Michele L. and David D. Pollard, 1997. Bedding plane slip in initial stages of fault-related folding Journal of Structural Geology Special Issue on Fault-Related Folding, vol. 19, pp. 567-581.

Cooke, Michele L., 1997.  Fracture localization along faults with spatially varying friction, Journal of Geophysical Research, vol. 102, pp. 22,425-22,434.

Savage, Heather M. and Michele L. Cooke, 2010. Unlocking the effects of friction on fault damage zones, Journal of Structural Geology, 1732-1741, Description: cleardoi:10.1016/j.jsg.2009.08.014.

http://www.sciencedirect.com/science/article/pii/S0191814114000959

McBeck, Jessica, Elizabeth H. Madden and Michele L. Cooke, 2016. Growth by Optimization of Work (GROW): A new modeling tool that predicts fault growth through work minimization, Computers and Geoscience, 10.1016/j.cageo.2015.12.019

## Authors

**Kelian Dascher-Cousineau** workflow, **Michele Cooke** - Fric2d, **Jess Mecbeth** GROW

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
