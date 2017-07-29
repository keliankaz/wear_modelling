function [] = fric2D(fileName,X,Y)
% Automated General work flow for the fric2d software 

% INPUT:
% filename:     string specifying the name of the fric2D input file (no
%               extension)
% X:            N by 1 array of x-coordinates of a desired fault profile
% Y:            N by 1 array of y-coordinates of a desired fault profile
% 
% varargin:
% not yet set...

% OUTPUT:
% 

% must be in same folder as:
% - simple_shear_bc.pl
% - fric2D files
% this file create the boundary line geometry and boundary stress
% conditions (this might be changed in future edits so that all inputs
% are set within the matlab environment - specifically for rock properties)

% must have the following paths to programs
% to fric2d (version fric2d_3.2.8):         fric2d_code/fric2d_3.2.8 -i
% to GROW (version GROW_JULY_2016.pl):      ./GROW_JULY_2016.pl

% TO CHANGE...                       (for the moment)
% ------------------------------------------------------------------------

% boundary conditions:      
% simple_shear_bc.pl

% fault geometery: 
% X Y inputs

% padding(size of the box relative to the fault): 
% makeinputfile: 176 ->                     [meanSegmentLength,boxSize] = makeprofile(profileName,X,Y,2,40);
                            
% endbit displacement:
% makeinputfile: 209-210 ->                 BVSTop = -0.00002; BVSBot = -0.0002;

% failure parameters: cohesion and angle of internal frinction
% assessfailure: 458-459                    c           = 5;    % cohesion
%                                                         theta       = 20;   % angle of internal friction


%--------------------------------------------------------------------------
% system commands are for set for usage with mac or linux
% -------------------------------------------------------------------------

%% input parameters:
% -----------------------------------------------------------------------
inputParameters  = userInput();
% -----------------------------------------------------------------------

%% work flow:

% *check Input*:
[X, Y] = checkinput(fileName, X,Y,nargin);

%% 1) make input file
[meanSegmentLength, newX] = makeinputfile(fileName,X,Y,inputParameters);

%% 2) run program 
runfric2D(fileName);

%% 3) extract info
outputStruct = extractoutput(fileName);

%% 4) assess shear failure
[shearFailureCoord  , tensileFailureCoord, ...
 coulombCriterion   , tensileCriterion]   = assessfailure(newX,Y, outputStruct, inputParameters);
                               
%% 5) feed into GROW
    
% 5a)choose which crack(s) to grow (now set to 'all', could alternatively
% be one of: 'all', 'tensile', 'shear', 'max coulomb', and 'max tensile')
crackLocations = choosecracks(shearFailureCoord  , tensileFailureCoord  , ...
                              inputParameters.user.cracks2grow, inputParameters);

% 5b) run GROW with the chosen ckracks
growOutputFile = runGROW(fileName, crackLocations, meanSegmentLength, inputParameters);    

%% 6) extract GROW output (not finished)
% GROWOutputStruct = extractGROWoutput(growOutputFile, crackLocations);

%% 7)  graphically show entire output
plotoutputinfo

%% 8) organize files
% create directory where all files will be stored

organizefiles

% ------------------------------------------------------------------------------------------------------
%% embeded functions 

% 7) plot output information into three plots
    function [] = plotoutputinfo(varargin)
    
        if ~strcmp(inputParameters.user.showFric2DOutput,'yes')
            if strcmp(inputParameters.user.showFric2DOutput, 'no')
                disp('output not displayed change inputParameters.user.showFric2DOutput in input section to ''yes'' if so desired')
            else
                error('inputParameters.user.showFric2DOutput in input section must be ''yes'' or ''no''')
            end
            
            return
            
        end
        
% three plots: 

% a) geometry
% b) displacements
% c) stresses
% d) failure criterion

%%%%%%%%%%%%%%%%%%%%%%%%%

fric2dFig = figure;
numPlots        = 4;
subPlotCount    = 1;

% a) geometry: (boundary box, fault, failure points)
subplot(numPlots, 1, subPlotCount) 
subPlotCount = subPlotCount+1;
hold on

boundaryBox     = fill( outputStruct.boundaryLines.xBeg , ...
                                   outputStruct.boundaryLines.yBeg , ...
                        [0.8 0.8 0.8]);
faultProfile                = plot( outputStruct.faultInfo.xBeg     , ...
                                             outputStruct.faultInfo.yBeg     , ...
                                             '.-');
leftEndBitLine           = plot( outputStruct.leftEndBit.xBeg, ...
                                             outputStruct.leftEndBit.yBeg, 'r.-');
rightEndBitLine         = plot( outputStruct.rightEndBit.xBeg, ...
                                             outputStruct.rightEndBit.yBeg, 'r.-');
                   
tensileFailurePoints    = scatter(tensileFailureCoord(:,1), tensileFailureCoord(:,2));
shearFailurePoints      = scatter(shearFailureCoord(:,1), shearFailureCoord(:,2));

hLegendA        = legend([boundaryBox           , ...
                          faultProfile          , ...
                          rightEndBitLine        , ...
                          tensileFailurePoints  , ...
                          shearFailurePoints   ], ...
                          ...
                          'Boundary Elements'   , ...
                          'Rough Fault Profile' , ...
                          'end bits'            , ...
                          'Failure in tension'  , ...
                          'Failure in shear'    ); 
                     
set( gca                                        , ...
     'Box'          , 'off'                     , ...
     'FontName'     , 'Helvetica'               );
axis equal
     
% b) displacmement (shear, normal)
subplot(numPlots, 1, subPlotCount) 
subPlotCount = subPlotCount+1;
hold on  

shearDisplacementPlot   = plot(outputStruct.faultInfo.xBeg, ...
                               outputStruct.faultStress.DS);                     
normalDisplacementPlot  = plot(outputStruct.faultInfo.xBeg, ...
                               outputStruct.faultStress.DN);
                           
hLegendB        = legend([shearDisplacementPlot , ...
                          normalDisplacementPlot],...
                          ...
                          'Shear displacement'  , ...
                          'Normal displaceement');
                      
% c) stresses (shear, normal, tangential)
subplot(numPlots, 1, subPlotCount) 
subPlotCount = subPlotCount+1;
hold on

shearStressPlot         = plot(outputStruct.faultInfo.xBeg, ...
                               outputStruct.faultStress.SigmaS);
normalStressPlot        = plot(outputStruct.faultInfo.xBeg, ...
                               outputStruct.faultStress.SigmaN);
                           
tangentialStressPlot    = plot(outputStruct.faultInfo.xBeg(2:end-1), ...
                               outputStruct.faultStress2.Tangential);
                           
hLegendC        = legend([shearStressPlot       , ...
                          normalStressPlot      , ...
                          tangentialStressPlot],...
                          ...
                          'Shear Stress'        , ...
                          'Normal Stress'       , ...
                          'Tangential Stress'   );

% d) failure criterion (coulomb, tensile)
subplot(numPlots, 1, subPlotCount) 
subPlotCount = subPlotCount+1;
hold on

coulombPlot             = plot(outputStruct.faultInfo.xBeg(2:end-1), ...
                               coulombCriterion);
tensileStressPlot       = plot(outputStruct.faultInfo.xBeg(2:end-1), ...
                               tensileCriterion);
hLegendD        = legend([coulombPlot           , ...
                          tensileStressPlot     ], ...
                          ...
                          'coulomb failure criterion (\tau_m - (\sigma_m sin\phi+c cos\phi) > 0)', ...
                          'tensile failure criterion (\sigma_3 - c > 0)');
    end

% 8) organize files to reduce clutter in directory
    function organizefiles
    
% make new directory for the files from this run
newDirectoryName    = ['output/', fileName];

if ~isdir(newDirectoryName)
command             = sprintf('mkdir %s', newDirectoryName);
system(command);
end 

% move all file that start with 'fileName' into new directory
command = sprintf('mv %s %s', [fileName,'.*'], [newDirectoryName,'/']);
system(command);
command = sprintf('mv %s %s', 'grow_specs.*', [newDirectoryName,'/']);
system(command);

% create input and output directories

if ~isdir([newDirectoryName, '/inputs'])
command = sprintf('mkdir %s', [newDirectoryName, '/inputs']);
system(command);
end

if ~isdir([newDirectoryName, '/outputs'])
command = sprintf('mkdir %s', [newDirectoryName, '/outputs']);
system(command);
end

if ~isdir([newDirectoryName, '/efficient_files'])
command = sprintf('mkdir %s', [newDirectoryName, '/efficient_files']);
system(command);
end

if ~isdir([newDirectoryName, '/prev'])
command = sprintf('mkdir %s', [newDirectoryName, '/prev']);
system(command);
end

% move inpputs an outputs to their respective locations
command = sprintf('mv %s %s', [newDirectoryName,'/*.in'], [newDirectoryName,'/inputs']);
system(command);
command = sprintf('mv %s %s', [newDirectoryName,'/*.out'], [newDirectoryName,'/outputs']);
system(command);
command = sprintf('mv %s %s', [newDirectoryName,'/*.eff'], [newDirectoryName,'/efficient_files']);
system(command);
command = sprintf('mv %s %s', [newDirectoryName,'/*.prev'], [newDirectoryName,'/prev']);
system(command);

profileName                 = [fileName,'_profile.txt'];
command = sprintf('mv %s %s', [newDirectoryName,'/*.prev'], [newDirectoryName,'/prev']);
system(command);

system('rm end_bits.txt');
system('rm input');
profileName = [fileName,'_profile.txt']; system(['rm ',profileName]);

    end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------- User Input --------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inputParameters = userInput()
    inputParameters = [];

    % material parameters:

    inputParameters.material.E                      = 30000;     % Young's Modulus MPa
    inputParameters.material.nu                    = 0.25;      % Poisson's ratio

    inputParameters.material.c                      = 5;         % cohesion (MPa)
    inputParameters.material.phi                   = 20;        % angle of intenral friction  (degrees)

    % boundary conditions

    inputParameters.boundaryConditions.sigxx       = 100;       % normal stress parallel to faults MPa
    inputParameters.boundaryConditions.sigyy        = 100;       % normal stress perpendicular to faults MPa
    inputParameters.boundaryConditions.shear       = 40;        % shear stress on faults  MPa

    % end bits
    inputParameters.boundaryConditions.BVSTop       = -0.0000445;  % displacement on left end bit (m)
    inputParameters.boundaryConditions.BVSBot       = -0.000215;   % displacement on right end bit (m) 

    % user choices

    inputParameters.user.showFric2DOutput           = 'no';     % plot output from fric2d ('yes' or 'no')
    inputParameters.user.cracks2grow                = 'max coulomb';     % 'all':           all the coordinates specified in input (tensile+shear)
                                                                 % 'tensile':       only cracks that failed in tensile stress (listed in
                                                                 %                  tensileFailureCoord)
                                                                 % 'shear'          only cracks that failed in shear stress (listed in
                                                                 %                  shearFailureCoord)
                                                                 % 'max tensile'    most tensile crack
                                                                 % 'max coulomb'      most sheared crack                                                            
    inputParameters.user.runGROWModel               = 'yes' ;     % run grow ('yes' or 'no')
    % GROW inputs are all in the runGROW function
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------- section functions ----------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input
function [X,Y] = checkinput(fileName, X,Y, numInput)
        
        if isa(fileName, 'string')
            error('input '' fileName must be of ype ''string''');
        end
    
        if numInput ~= 3
            error('There must be three inputs: an X array, a Y array and a string specifying the file name')
        end
        
        if min(size(X)) ~=1
            error('X array must be 1 by N or N by 1')
        end
        
        if min(size(Y)) ~=1
            error('Y array must be 1 by N or N by 1')
        end
        
        if length(X) ~= length(Y)
            error('X and Y arrays must be same length')
        end
        
        if length(X(1,:)) == length(X)
            X = X';
        end
        
         if length(Y(1,:)) == length(Y)
            Y = Y';
         end
         
    end

%% 1) make the input file for the fric2D program
function [meanSegmentLength, X] = makeinputfile(fileName,X,Y, inputParameters)

% create profile: 
profileName                 = [fileName,'_profile.txt']; % dont change this same name used in 6b

padding = 2;
boxAspectRatio = 0.2;

[X, meanSegmentLength,boxSize] = makeprofile(profileName,X,Y,padding, boxAspectRatio);

% create boundary condition
    % input ptspacing based on number of points
    % simple_shear_bc.pl 0.06 0.15 0.0002
    % (based on profile1)
    
commandFormat = 'perl simple_shear_bc.pl %f %f %f %f %f %f %f %f';    
command = sprintf(commandFormat         , ...
            boxSize(2)                  , ...
            boxSize(1)                  , ...
            meanSegmentLength           , ...
            inputParameters.material.E  , ...
            inputParameters.material.nu , ...
            inputParameters.boundaryConditions.sigxx, ...
            inputParameters.boundaryConditions.sigyy, ...
            inputParameters.boundaryConditions.shear);
            
system(command);

% make the endbit so that there is not a displacement discontinuity at the
% end of the fault 
BVSTop      = inputParameters.boundaryConditions.BVSTop;
BVSBot      = inputParameters.boundaryConditions.BVSBot;

% calculate number of endbit segments on either side of input fault
% geometry (magic number 4 is used to calulate half the distance from the
% tip of the fault geometry to the end  side of the box)
numBufferSegments = floor(boxSize(1)/(padding(1)*4*meanSegmentLength));
make_end_bits(X,Y, numBufferSegments, BVSTop, BVSBot);
% make the length of the end bit a function of the profile length

% concatenatem to the input file
command = sprintf('cat %s %s %s > %s', 'input', profileName, 'end_bits.txt', [fileName, '.in']);
system(command);

% more options...
    % options that may need to be tweeked: 
    % observation points (change 'odata' and 'nolines')
    % tolerance
    
% 1a) make the fault geometry based on series of x-y points, get the mean BEM
% size and determine the appropriate box size
function [newX, MEAN_SEGMENT_LENGTH,boxSize] = makeprofile(fileName,X,Y,padding, boxAspectRatio)

% prints profile defined by vertical vectors X and Y in line element form
% in the file specified by input fileName
% padding [dim1,dim2] defines the size of the boundary condition box as a function of
% the profile length and amplitude of the profile.

% *edits*:
% wed jul12 2017 changed the padding so that vertical box dimensions are
% instead a function of the profile length and determined by  the
% 'boxAspectRatio'

RX = range(X);
RY = range(Y);

MEAN_SEGMENT_LENGTH = mean(sqrt((X(2:end)-X(1:(end-1))).^2 +(Y(2:end)-Y(1:(end-1))).^2));

% boxSize = padding.*[RX,RY]; % not great because the box size changes as a
% function of topography, which is bad if set in a loop

% instead:

boxSize = zeros(1,2); 
boxSize(1) = padding(1)*RX;
boxSize(2) = boxSize(1)*boxAspectRatio;

% errors/warnings related to profile in relation to box size
if max(abs(Y))>boxSize(2)*0.1
    if max(abs(Y))>boxSize(2)/2
        error('Profile topography is too large and intersect the boundary box, choose a larger aspect ratio for the box by changing the vatiable ''boxAspectRatio'' ')
    end
    warning('Profile topography is large, consider choosing a larger aspect ratio for the box by changing the vatiable ''boxAspectRatio'' ');
end

% roubnd to increments of the length of elements

newX = X - min(X) + boxSize(1)/2 - RX/2; % shift X to be centered in the observation box

fid=fopen(fileName,'w');

for n = 1:(length(newX)-1)
    fprintf(fid, ...
        '1 %f %f %f %f 0 1.00E+10 0 0 0 1 1 0\n', ...
        [newX(n),  Y(n),  newX(n+1),   Y(n+1)]');
end

fprintf(fid, '\n');

fclose(fid);

end

% 1b) add end bits to input file according to numBuffer segments (figiting
% needed)
function [] = make_end_bits(X,Y,numBufferSegments, BVSTop, BVSBot)

    intro =    {'*FRACTURES'; ...
    '*num	xbeg	ybeg	xend	yend	kode	bvs     bvn     bvs     bvn'; ...
    '*----	----	----	----	----	----	------	----	-----	--------'};

fileID = fopen('end_bits.txt','w');
fprintf(fileID,'%s\r\n',intro{:});

X1      = X(1);
Y1      = Y(1);
Xend    = X(end);
Yend    = Y(end);
start1  = [X1-numBufferSegments*meanSegmentLength,Y1];
end1    = [X1,Y1];
start2  = [Xend,Yend];
end2    = [Xend+numBufferSegments*meanSegmentLength,Yend];

% put in a formulation so that the velocity boundaries are defined be a
% strain rate and the length and height of the box rather than arbitrary
% numbers
% strainRate =   ....

%End bits:      left end bit    right end bit
%               -------------   ---------------
kode        = [ 3,              3               ];
bvs         = [ 0,              0               ];
bvn         = [ 0,              0               ];
BVS         = [ BVSTop,         BVSBot          ];
BVN         = [ -100;           -100            ];

% make the rows of data for input file

ROW1    = [numBufferSegments, start1,end1,kode(1),bvs(1),bvn(1),BVS(1),BVN(1)];
ROW2    = [numBufferSegments, start2,end2,kode(2),bvs(2),bvn(2),BVS(2),BVN(2)];

% write them to the file
fprintf(fileID, 'fracture left no no');
fprintf(fileID, '\n %f %f %f %f %f %f %f %f %f %f \n\n',ROW1);
fprintf(fileID, 'fracture right no no');
fprintf(fileID, '\n %f %f %f %f %f %f %f %f %f %f \n\n',ROW2);
fclose(fileID);

end

end

%% 2) run fric2D
function runfric2D(fileName)
% command = ['fric2d_3.2.8 -i ', ...
%            fileName,'.in -v -o ', fileName,'.out'];
% system(command);

command = ['./fric2d -i ', ...
           fileName,'.in -v -o ', fileName,'.out'];
system(command);

end

%% 3) extract output into a structure (uses getsection, assign2struct)
function [OUTPUT_STRUCT] = extractoutput(fileName)

% extracts numeric blocks of information out of output file and saves them
% into matlab structure with fields names corresponding (sometimes with
% small adjustments to header lines)

% input:
% fileName: name of the file in a string without the extension (here
% assumed to be  .out). IMPORTANT: output file is assumed to have 9 data
% blocks. Namely 2 boudaryu line blocks, followed by 4 'end bit' blocks and
% 3 fault info blocks. The current state does not extract 'end bits'
% information.

% output:
% structure with with fields corresponding to blocks of data:
% 1) boundaryLines
% 2) boundaryConditions
% 3) faultinfo:     information about fault geometry etc.
% 4) faultStress
% 5) faultStress2:  additional fault stresses namely including the
%                   tangential stresses

% function calls:
% read_blocks_of_numeric_data
% assign2struct

% usage:
% extractoutput(sample_section_test)

% 2017-04-17: tested and functions properly using the file
% sample_section_test.out


outputFileName = [fileName,'.out'];

% load in block of data (TERRIBLE IDEA BUT FUCK IT): 
expectedDataBlocks = 9; % this is here to make sure that I'm not fucking up the parsing of the output

dataBlocks = read_blocks_of_numerical_data(outputFileName, 1000);
numBlocks  = length(dataBlocks);

if numBlocks ~= expectedDataBlocks
    error('the output is not of the expected size consider double-checking output parsing')
end

%extract the boundary lines:
boundaryLineBlock       =     dataBlocks(1);
boundaryLineHeaders     =     {'Elt'            , ...
                               'xBeg'           , ...
                               'yBeg'           , ...
                               'xEnd'           , ...
                               'yEnd'           , ...
                               'length'         , ...
                               'angle'          , ...
                               'USstat'         , ...
                               'UNstat'         , ...
                               'USMonotonic'    , ...
                               'UNMonotonic'    };

boundaryLines           =       assign2struct(boundaryLineBlock{:},boundaryLineHeaders);                      

% extract boundary conditions:
boundaryConditionBlock  = 	dataBlocks(2);
boundaryConditionHeaders = {    'Elt'           , ...
                                'DS'            , ...
                                'USminus'         , ...
                                'USplus'         , ...
                                'DN'            , ...
                                'UNminus'         , ...
                                'UNplus'         , ...
                                'SigmaS'        , ...
                                'SigmaN'        };

boundaryConditions =        assign2struct(boundaryConditionBlock{:},boundaryConditionHeaders);

% exctract fault information

faultInfoBlock          =       dataBlocks(7);
faultInfoHeaders        =       boundaryLineHeaders;

faultInfo               =       assign2struct(faultInfoBlock{:},    faultInfoHeaders);

% exctract fault information: stress conditions

faultStressBlock        =       dataBlocks(8);
faultStressHeaders      =       boundaryConditionHeaders;

faultStress             =       assign2struct(faultStressBlock{:},  faultStressHeaders);

% exctract fault information: stress conditions

faultStress2Block       =       dataBlocks(9);
faultStress2Headers     =    {  'BE'           , ...
                                'Tangential'    , ...
                                'maxPrinctop'   , ...
                                'angleTop'      , ...
                                'tangential'    , ...
                                'maxPrincBot'   , ...
                                'angleBot'      };

faultStress2            =     assign2struct(faultStress2Block{:},   faultStress2Headers);


% extract endbit geometry:
leftEndBitBlock         =       dataBlocks(3);
rightEndBitBlock        =       dataBlocks(5);
endBitHeaders           =       boundaryLineHeaders;

leftEndBit              =       assign2struct(leftEndBitBlock{:},   endBitHeaders);
rightEndBit             =       assign2struct(rightEndBitBlock{:},  endBitHeaders);


OUTPUT_STRUCT                   = [];
OUTPUT_STRUCT.boundaryLines     = boundaryLines;
OUTPUT_STRUCT.boundaryConditions= boundaryConditions;
OUTPUT_STRUCT.faultInfo         = faultInfo;
OUTPUT_STRUCT.faultStress       = faultStress;
OUTPUT_STRUCT.faultStress2      = faultStress2;
OUTPUT_STRUCT.leftEndBit        = leftEndBit;
OUTPUT_STRUCT.rightEndBit       = rightEndBit;



end

%% 4) assess failure based on model output
function [shearFailureCoord,tensileFailureCoord, ...
          coulombCriterion, tensileCriterion] = assessfailure(X, Y, OUTPUT_DATA_STRUCT, ...
                                                              inputParameters)

% assessed the failure of BE's along the fault profile. This is done by 1)
% using the structure output produced by the extractouput function to
% determine the stress tensor, 2) finding the prinple stesses (the
% eigenvalues of the stress tensor) and 3) assessing both shear
% (mohr-coulomb) and tensile failure criterions. 4) finding the
% corresponding locations on the profile

% INPUT:

% OUTPUT_DATA_STRUCT: output data extracted using extractoutput

% OUTPUT:
% shearFailureCoord:    2 by N array with x-y coordinates of elements that failed in shear        
% tensileFailureCoord:  2 by N array with x-y coordinates of elements that failed in tnesion
% coulombCriterion:     1 by N array coulomb stress 
%                       (tau_m-[\sigma_m*sin(\phi) + c*cos(\phi)])
%                       where
%                       \sigma_m = (\sigma_1+\sigma_3)/2
%                       \tau_m = (\sigma_1-\sigma_3)/2
% tensileCriterion:     1 by N array tensile criterion
%                       (-\sigma_3-c)

% if the latter 2 are positive, the fault elements are expected to fail

% function calls:
% gettensor

% usage:
% assessfailure(outputDataStruct)

% 2017-04-17: tested and functions properly using the file
% sample_section_test.out processed with function extractoutput (did not
% test the physical viability of results)


c           = inputParameters.material.c;    % cohesion 
phi         = inputParameters.material.phi;   % angle of internal friction
      
[sigma11, sigma22, sigma12] = gettensor(OUTPUT_DATA_STRUCT,'top');
[shearFailureCoordTop,    tensileFailureCoordTop]     = find_failure_points;

[sigma11, sigma22, sigma12] = gettensor(OUTPUT_DATA_STRUCT,'bottom');
[shearFailureCoordBottom, tensileFailureCoordBottom]  = find_failure_points;

shearFailureCoord   = [shearFailureCoordTop; shearFailureCoordBottom];
tensileFailureCoord = [tensileFailureCoordTop; tensileFailureCoordBottom];

function [shearFailureCoord, tensileFailureCoord] = find_failure_points()

% returns sorted (as a function of magnitude) coordinates for failure
% both in shear and in tensile stress) 

% uses variables: sigma11, sigma22, sigma12, OUTPUT_DATA_STRUCT
    
numBE           = length(sigma11);

coulombCriterion= zeros(numBE,1);
tensileCriterion= zeros(numBE,1);

for iBE = 1:numBE
    
    tensor      =  [sigma11(iBE), sigma12(iBE); ...
                    sigma12(iBE), sigma22(iBE)];
    eigenValues = eig(tensor);
    sigma1      = eigenValues(1);
    sigma3      = eigenValues(2);
    
    tau_m       = (sigma1-sigma3)/2;
    sigma_m     = (sigma1+sigma3)/2;
    
    coulombCriterion(iBE) = tau_m - (sigma_m*sind(phi)+c*cosd(phi));
    tensileCriterion(iBE) = -sigma3-c;
 
end

% find points of failure (findpeaks is used to find local maxima in case
% there are segments of points that are above failure)
tempCC = coulombCriterion;
tempTC = tensileCriterion;

tempCC(tempCC<0) = 0; % remove section that are not going to fail
tempTC(tempTC<0) = 0;

[~,shearFailureInd  ] = findpeaks(tempCC,'SortStr', 'descend'); % find peak stresses and sort their indices
[~,tensileFailureInd] = findpeaks(tempTC,'SortStr', 'descend');

shearFailureCoord     =  [X(shearFailureInd+1), ... 
                                        Y(shearFailureInd+1)];
                               
tensileFailureCoord   =  [X(tensileFailureInd+1), ...
                                       Y(tensileFailureInd+1)];

                     
                     
end

end

%% 5a) choose which cracks to grow
function   CRACKLOCATION = choosecracks(shearFailureCoord  , tensileFailureCoord  , ...
                                        userChoice, inputParameters)
 % shoose which crack to grow:
 % 'all':           all the coordinates specified in input (tensile+shear)
 % 'tensile':       only cracks that failed in tensile stress (listed in
 %                  tensileFailureCoord)
 % 'shear'          only cracks that failed in shear stress (listed in
 %                  shearFailureCoord)
 % 'max tensile'    most tensile crack
 % 'max shear'      most sheared crack
 
% check if this function should be run
if ~strcmp(inputParameters.user.runGROWModel, 'yes')
    
    % if set to inputParameters.user.runGROWModel is set to 'no', indicate
    % that grow model is not lauched and where this can be changed
    if strcmp(inputParameters.user.runGROWModel, 'no')
    disp('No GROW Model was lauched, change ''inputParameters.user.runGROWModel'' input to ''yes'' if so desired')
    
    % send error is inputParameters.user.runGROWModel is not set to 'yes'
    % or 'no'
    else;  error('inputParameters.user.runGROWModel in input section must be ''yes'' or ''no''')
    end
    
    % leave function if inputParameters.user.runGROWModel is not 'yes'
    return
end


if      strcmp(userChoice, 'all')
    CRACKLOCATION = [shearFailureCoord;tensileFailureCoord];
    warning('Running all crack locations may be computationally expensive')
elseif  strcmp(userChoice, 'tensile')
    CRACKLOCATION = tensileFailureCoord;
    
elseif  strcmp(userChoice, 'shear')
    CRACKLOCATION = shearFailureCoord;
    
elseif  strcmp(userChoice, 'max tensile')
    CRACKLOCATION = tensileFailureCoord(1,:); % max is the first location because the list is sorted
    
elseif  strcmp(userChoice, 'max coulomb')
    CRACKLOCATION = shearFailureCoord(1,:); % max is the first location because the list is sorted
else
    error('choosecracks ''userChoice'' input must be one of: ''all'', ''tensile'',''shear'',''max tensile'' or ''max coulomb''')
end
                  
end

%% 5b) use outputs from fric2D and failure asseessment to track fracture growth
function growFileName = runGROW(fileName,crackLocations, BELength, inputParameters)
% this function runs GROW (GROwth by Work-minimization)

% input:
% fileName: string with the name of a fric2D input file 

% check if the function was called
if ~strcmp(inputParameters.user.runGROWModel, 'yes')
    if strcmp(inputParameters.user.runGROWModel, 'no')
        disp('No GROW Model was lauched, change ''inputParameters.user.runGROWModel'' input to ''yes'' if so desired')
    else 
        error('inputParameters.user.runGROWModel in input section must be ''yes'' or ''no''')
    end
    return
end

% grow input file
growFileName = 'grow_specs.in';
path2GrowFile = 'GROW/';

% grow input  parameters parameter

% this will be part of the command to run GROW
angleResolution         = 45;
startAngle              = 45;
endAngle                = 315;

% parameters tacked on to input file
fault_flaw              = [];
fault_flaw.tag          = 'Flaw-Fault'; % tag
fault_flaw.stiffS       = 10^10;        % shear stiffness           
fault_flaw.stiffN       = 10^10;        % normal stiffness    
fault_flaw.cohes        = 0;            % sliding cohesion
fault_flaw.friction_s   = 0.6;          % static friction
fault_flaw.friction_d   = 0.0;          % dynamic friction
fault_flaw.L            = 10^-5;        % critical sliding distance


intact_flaw             = [];
intact_flaw.Tag         = 'Flaw-intact';
intact_flaw.fault       = 'crittical_point';
intact_flaw.xcoord      = crackLocations(:,1);
intact_flaw.ycoord      = crackLocations(:,2);
intact_flaw.length      = BELength;             
intact_flaw.grow_both   = 'yes';
intact_flaw.stiffS      = 0;
intact_flaw.stiffN      = 10^10;
intact_flaw.T           = 5;
intact_flaw.S_0         = 50;
intact_flaw.cohesion    = 0;
intact_flaw.friction_s  = 0.6;
intact_flaw.friction_d  = 0.0;
intact_flaw.L           = 10^-5;

tempFile = 'grow_temp';
fid=fopen(tempFile,'w');

for iFlaw = 1:length(intact_flaw.xcoord)

% growInputSpecL1 = '*%s %s %s %s %s %s %s \n';
% fprintf(fid,growInputSpecL1, ...
%                                 'tag', ...
%                                 'stiffS', ...
%                                 'stiffN',...
%                                 'cohes', ...
%                                 'fric-s', ...
%                                 'fric-d',...
%                                 'L');
% 
% growInputSpecL2 = '*%s %f %f %f %f %f %f \n';
% fprintf(fid,growInputSpecL2, ...
%                                 fault_flaw.tag              , ...
%                                 fault_flaw.stiffS           , ...
%                                 fault_flaw.stiffN           , ...
%                                 fault_flaw.cohes            , ...
%                                 fault_flaw.friction_s       , ...
%                                 fault_flaw.friction_d       , ...
%                                 fault_flaw.L                );
%                             
% growInputSpecL3 = '*%s %s %s %s %s %s %s %s %s %s %s %s %s %s \n';                            
% fprintf(fid,growInputSpecL3, ...
%                                 'Tag', ...
%                                 'fault', ...
%                                 'xcoor', ...
%                                 'ycoor', ...
%                                 'length', ...
%                                 'grow_both?', ...
%                                 'stiffS', ...
%                                 'stiffN', ...
%                                 'T', ...
%                                 'S_0', ...
%                                 'cohes', ...
%                                 'friction-s', ...
%                                 'friction-d', ...
%                                 'L');
% 
% growInputSpecL4 = '%s %s %f %f %f %s %f %f %f %f %f %f %f %f\n';
% fprintf(fid,growInputSpecL4, ...
%                                 intact_flaw.Tag             , ...
%                                 [intact_flaw.fault,num2str(iFlaw)], ...
%                                 intact_flaw.xcoord(iFlaw)   , ...
%                                 intact_flaw.ycoord(iFlaw)   , ...
%                                 intact_flaw.length          , ...
%                                 intact_flaw.grow_both       , ...
%                                 intact_flaw.stiffS          , ...
%                                 intact_flaw.stiffN          , ...
%                                 intact_flaw.T               , ...
%                                 intact_flaw.S_0             , ...
%                                 intact_flaw.cohesion        , ...    
%                                 intact_flaw.friction_s      , ...    
%                                 intact_flaw.friction_d      , ...
%                                 intact_flaw.L               );

growInputSpecL1 = '*%s %s %s %s %s %s %s %s %s %s %s %s %s %s \n';                            
fprintf(fid,growInputSpecL1, ...
                                'Tag', ...
                                'fault', ...
                                'xcoor', ...
                                'ycoor', ...
                                'length', ...
                                'grow_both?', ...
                                'stiffS', ...
                                'stiffN', ...
                                'T', ...
                                'S_0', ...
                                'cohes', ...
                                'friction-s', ...
                                'friction-d', ...
                                'L');

growInputSpecL2 = '*%s %f %f %f %f %f %f \n';
fprintf(fid,growInputSpecL2, ...
                                fault_flaw.tag              , ...
                                fault_flaw.stiffS           , ...
                                fault_flaw.stiffN           , ...
                                fault_flaw.cohes            , ...
                                fault_flaw.friction_s       , ...
                                fault_flaw.friction_d       , ...
                                fault_flaw.L                );
                           

growInputSpecL4 = '*%s %s %f %f %f %s %f %f %f %f %f %f %f %f\n\n';
fprintf(fid,growInputSpecL4, ...
                                intact_flaw.Tag             , ...
                                [intact_flaw.fault,num2str(iFlaw)], ...
                                intact_flaw.xcoord(iFlaw)   , ...
                                intact_flaw.ycoord(iFlaw)   , ...
                                intact_flaw.length          , ...
                                intact_flaw.grow_both       , ...
                                intact_flaw.stiffS          , ...
                                intact_flaw.stiffN          , ...
                                intact_flaw.T               , ...
                                intact_flaw.S_0             , ...
                                intact_flaw.cohesion        , ...    
                                intact_flaw.friction_s      , ...    
                                intact_flaw.friction_d      , ...
                                intact_flaw.L               );

end

fclose(fid);

profileName                 = [fileName,'_profile.txt'];
command = sprintf('cat %s %s %s %s > %s', 'input', profileName, tempFile, 'end_bits.txt', [growFileName]);
system(command);

command = sprintf('rm %s', tempFile);
system(command);

command = sprintf('perl GROW_nov_17_15.pl %s %f %f %f'         , ...
                                growFileName                , ...
                                angleResolution             , ...
                                startAngle                      , ...
                                endAngle                    );
                            
system(command);


% stuff from Michele
%- add lines in the Fault conditions section: (make sure parameters
% are rigth for sandstone)

% *Tag fault xcoor ycoor length grow_both? stiffS stiffN T S_0 cohes friction-s friction-d L
% *Flaw-Fault 			 				   	   0	1e10		0	0	0	1	     	1			0.0	
% *Flaw-Intact 	rough	0.055424 0.01173379	0.0006	yes         1e10	  1e10		5	50	0	0.6	     	0			1e-5

% to run: ./GROW_JULY_2016.pl bump.in 45 45 315


end

%% 6) extract Grow output (being  built)

%% 7) visually diaplay output information (function is located in embeded function section)

%% 8) Organize files (function is located in embeded function section)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------- Being built ---------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% extract Grow output

function GROWOutputStruct = extractGROWoutput(crackLoc, growFileName)

% extracts numeric blocks of information out of GROW output file and saves
% them into matlab structure with fields names corresponding (sometimes
% with small adjustments to header lines)

% % input: % fileName: name of the file in a string without the extension
% (here % assumed to be  .out). IMPORTANT: output file is assumed to have 9
% data % blocks + x blocks per flaw. Namely 2 boundary line blocks,
% followed by 4 'end bit' blocks and % 3 fault info blocks. The current
% state does not extract 'end bits' % information.

% output:
% structure with with fields corresponding to blocks of data:

% function calls:
% read_blocks_of_numeric_data
% assign2struct

% usage:
% extractoutput('sample_section_test')

% 2017-04-17: tested and functions properly using the file
% sample_section_test.out


% check if this function should be run:
if strcmp(inputParameters.user.runGROWModel, 'yes')
    if strcmp(inputParameters.user.runGROWModel, 'no')
        disp('No GROW Model was lauched, change ''inputParameters.user.runGROWModel'' input to ''yes'' if so desired')
    else 
        error('inputParameters.user.runGROWModel in input section must be ''yes'' or ''no''')
    end
    return
end
        
outputFileName = [fileName,'.out'];

% load in block of data (TERRIBLE IDEA BUT FUCK IT): 

% determine howmany cracks where grown
crackLocDim = size(crackLoc);
numCracks = crackLocDim(1);
blocksPerCrack = 2; % make sure this is correct
expectedFric2dDataBlocks = 9; % this is here to make sure that I'm not fucking up the parsing of the output
expectedDataBlocks = expectedFric2dDataBlocks + blocksPerCrack*numCrack;

dataBlocks = read_blocks_of_numerical_data(outputFileName, 1000);
numBlocks  = length(dataBlocks);

if numBlocks ~= expectedDataBlocks
    error('the GROW output is not of the expected size consider double-checking output parsing')
end

% Get crack info:

for iCrack = 1:numCrack
   % exctract data for each crack
   
   % steps:
% define headers (what will they be?)

% get block of data
% get the geometry

% What I need to do here:
% figure out how many blocks of data will be produced from GROW

% for now the key blocks of data that I want are the crack geometry
% need to know what the header sections will be as well
end % for iCrack = 1:numCrack

end % function GROWOutputStruct = extractGROWoutput(crackLoc, growFileName)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------- handy tools -------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% extract block of numerical data in text file
function    out = read_blocks_of_numerical_data( filespec, block_size, delimiter )

%   Within a block all rows must have the same number of "columns". Test with str2num. 

%   2014-06-09, poi: Cannot handle the sample below, since the first line "1" is 
%               included in the block. 
%       1
%   255  255  255  255 
%   255  255  255  255 
%   255  255  255  255 
%   255  255  255  255 
%   255  255  255   99 
%   255  255   97   99 

%   block_size  (characters)
%   row_size    (numerical items)
    
%{
g = read_blocks_of_numerical_data( 'blocks_of_numerical_data.txt', 64 )
g = read_blocks_of_numerical_data( 'set1.txt', 64, '\t' )
g = read_blocks_of_numerical_data( 'project.txt', 8 )  
g = read_blocks_of_numerical_data( 'read_specific_lines_from_tex.txt', 64 )  
g = read_blocks_of_numerical_data( 'small_textfile.txt', 50 ); % fails
g = read_blocks_of_numerical_data( 'large_textfile.txt', 50 ); % fails
g = read_blocks_of_numerical_data( '140708_LO_03_140710112412.txt', 50 ); 
g = read_blocks_of_numerical_data( 'RajRaj.txt', 150 ); 
%}
    narginchk( 2, 3 )
    buffer  = fileread( filespec );
    
    if nargin == 2 
        del_xpr = '[ ]+';
        trl_xpr = '[ ]*';
    else
        del_xpr = ['([ ]*',delimiter,'[ ]*)'];
        trl_xpr = ['([ ]*',delimiter,'?[ ]*)'];
    end
    
    num_xpr = '([+-]?(\d+(\.\d*)?)|(\.\d+))';
    sen_xpr = '([EeDd](\+|-)\d{1,3})?';         % optional scientific E notation
    num_xpr = [ num_xpr, sen_xpr ];
    
    nl_xpr  = '((\r\n)|\n)';
    
    row_xpr = cat( 2, '(^|', nl_xpr, ')[ ]*(', num_xpr, del_xpr, ')*'   ...
                    , num_xpr, trl_xpr, '(?=',nl_xpr,'|$)'              ); 
    
    blk_xpr = ['(',row_xpr,')+'];
    
    blocks  = regexp( buffer, blk_xpr, 'match' );
    
    is_long = cellfun( @(str) length(str)>=block_size, blocks );
    
    blocks(not(is_long)) = [];
    
    out = cell( 1, length( blocks ) ); 
    for jj = 1 : length( blocks )
        out{jj} = str2num( blocks{jj} );        
    end
end

%% assign info from loaded file into a structure array
function [STRUCTUREOUPUT] = assign2struct(dataArray,cellArray)

structOut   = [];
for iColHeader = 1:length(cellArray)
    structOut.(cellArray{iColHeader}) = dataArray(:,iColHeader);
end

STRUCTUREOUPUT = structOut;

end

%% find the stress tensor information from the pre-processes output
function    [sigma11, sigma22, sigma12] = gettensor(outputStruct, section)
% extracts tensor info for each element from output struct of specified section

%   INPUT:
%   outputStruct:   is the strucutre that is created from the function extract
%                   output
%   section:        is a string specifying what section to get (as of now functions
%                   for 'top' or 'bottom'

%  OUTPUT:
%  TENSOR:          4 component 2d tensor

% extract shear and element-normal stresses:

sigmaN = outputStruct.faultStress.SigmaN;
sigmaS = outputStruct.faultStress.SigmaS;

% extract tangential normal stress
if strcmp(section, 'top')
    sigmaT = outputStruct.faultStress2.Tangential;
elseif strcmp(section, 'bottom')
    sigmaT = outputStruct.faultStress2.tangential;
else 
    error('input ''section'' must be either string ''top'' or ''bottom''')
end

% make stress values same length (sigmaN and sigmaS are 2 longer that the
% tangential stresses. Not too sure why this is the case)

sigmaN = sigmaN(2:end-1);
sigmaS = sigmaS(2:end-1);


sigma12 = sigmaS;

if sigmaN > sigmaT;     sigma11 = sigmaN;
                        sigma22 = sigmaT;  
else;                   sigma11 = sigmaT;
                        sigma22 = sigmaN;
end
end
% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------function no longer in use ------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% extract sections of text file using system commands
function getsection(fileName, anotherFileName, startString, endString)

% extract sections of text file using system commands

% input: 

% fileName:   input file name from which the section will be extracted 
% sartString: string that marks the beginning of the section to be selected
% endString:  string that marks the end of the section to be selected
% 
% note that the section should have clear column headers, and data section
%
% output:
%
% structure array with 3 fields: 
% 1) DESIRED_SECTION.data:      data section
% 3) 
% Usage:
% DESIRED_SECTION = getsection(fileName, startString, endString)

command1 = ['sed -n -e '',/'    , ...
           startString          , ...
           '/,/'                , ...
           endString            , ...
           'p'' '               , ...
           fileName             , ...
           ' > '                , ...
           anotherFileName      ];
       
         
system(command1)                        % grab section between startString and endString

end

%% grab all the data blocks from the output data: 
function [DATABLOCKS] = grabdata(fileName, numBlock, numColumnArray, numHeaderArray)

%(NOTE THAT THIS RELIES ON THE POUTPUT  STRUCTURE TO BEIDENTICAL EVERY TIME
%- I.E. NUMBER OF HEADER COLUMNS, NUMBER OF DATA BLOCK)

fid=fopen(fileName);
for i=1:numBlock
    fmt=repmat('%f',1,numColumnArray(i));  % build the format string for the number columns
    a(i)=textscan(fid,fmt,'headerlines',numHeaderArray(i),'collectoutput',1); % read section
end
fid=fclose(fid);

end

%% parse user input
function S = setVal(S, name, input)
% set user specified inputs to a default structure based on "pair-wise
% input"

% S         : is a structure to whom 'name' fields will be assigned if they exist in
%             the 'input'
% name      : cell array of all fields to possibly set
% input     : cell array of pair-wise input values

if length(input) >= 2
    
    specifier   = input(1:2:end);
    value       = input(2:2:end);
    
    for iName = 1:length(name)
        
        fieldName   = name(iName);
        ind         = strcmp(fieldName, specifier); % check every odd input
        fieldName   = char(fieldName);
        
        if any(ind)
            S.(fieldName) = value{ind}; %  assign field val to even input
        end
        
    end
end
end
% -------------------------------------------------------------------------