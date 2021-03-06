function GROWOutputStruct = fric2d_workflow(fileName,X,Y,varargin)
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
% GROWOUTPUT = fric2d_workflow(fileName,X,Y)  - (if one output)
% GROWOUTPUT: structure with with fields corresponding to blocks of data
%   segments: 1 by N (N is numuber of fault segments, or, equivalently, the
%   number of defects+1) structure with fields corresponding to .eff
%   headers for the crack information (num    xbeg    ybeg    xend    yend
%   stiffS stiffN ten_str init_coh slid_coh stat_fric dy_fric crit_slip)
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

%--------------------------------------------------------------------------
% system commands are for set for usage with mac or linux
% -------------------------------------------------------------------------

%% input parameters:
% -----------------------------------------------------------------------
[inputParameters, GROWInputParameters]  = userInput();
% -----------------------------------------------------------------------

%% work flow:

% *check Input*:
[X, Y,plotOption] = checkinput(fileName, X,Y,varargin, nargin);

% *check Directory*
check_directory(GROWInputParameters);

%% 1) make input file
[~, newX] = makeinputfile(fileName,X,Y,inputParameters);

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

if ~isempty(crackLocations)                           
% 5b) run GROW with the chosen ckracks
%runGROW(fileName, crackLocations, meanSegmentLength, inputParameters, GROWInputParameters);    
runGROW(newX, Y, crackLocations,inputParameters, GROWInputParameters);    

%% 6) extract GROW output (not finished)
GROWOutputStruct = extractGROWoutput(crackLocations, GROWInputParameters.growFileName, inputParameters);

end

%% 7)  graphically show entire output
plotoutputinfo;

%% 8) organize files
% create directory where all files will be stored

organizefiles;

%% 9) check model - issue potential warnings:

checkmodel(X,Y,inputParameters, outputStruct, GROWOutputStruct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------------- embeded functions ------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

figure

longPlotOn        = strcmp(plotOption,'long plot');
numPlots        = 4;
subPlotCount    = 1;

if longPlotOn
subplot(numPlots, 1, subPlotCount)
end
subPlotCount = subPlotCount+1;

% a) geometry: (boundary box, fault, failure points)
hold on

% plot the poundary box as a rectangle
boundaryBox     = fill(             outputStruct.boundaryLines.xBeg , ...
                                    outputStruct.boundaryLines.yBeg , ...
                                    [0.8 0.8 0.8]);
                                
% make quiver plot of boundary conditions
angle   = outputStruct.boundaryLines.angle;
US      = outputStruct.boundaryLines.USMonotonic;
UN      = outputStruct.boundaryLines.USMonotonic;

hold on

% make vectors
elementDir =[cosd(angle),sind(angle)]./repmat((sqrt(cosd(angle).^2+sind(angle).^2)),1,2); %unit vecto in the direction of the element;
USVec = elementDir.*repmat(US,1,2);
R = [0 -1; 1, 0];
UNVec = (R*elementDir'.*repmat(UN,1,2)')'; % rotated 90 degrees

USQuivers = quiver(outputStruct.boundaryLines.xBeg , outputStruct.boundaryLines.yBeg , ...
                   USVec(:,1), USVec(:,2), ...
                   'AutoScale', 'off');

UNQuivers1 = quiver(outputStruct.boundaryLines.xBeg , outputStruct.boundaryLines.yBeg , ...
                   UNVec(:,1), UNVec(:,2), ...
                   'AutoScale', 'off');
               
UNQuivers2 = quiver(outputStruct.boundaryLines.xBeg , outputStruct.boundaryLines.yBeg , ...
                   -UNVec(:,1), -UNVec(:,2), ...
                   'AutoScale', 'off');
                                
% plot the fault segments outputed from grow (if clause is within the function)                    
plotGROWOutput('geometry')

% superpose the fault geometry
faultProfile                = plot( outputStruct.faultInfo.xBeg     , ...
                                   	outputStruct.faultInfo.yBeg     , ...
                                             '.-');
                                         
% plot endbits                                         
leftEndBitLine          = plot(     outputStruct.leftEndBit.xBeg, ...
                                    outputStruct.leftEndBit.yBeg, 'r.-');
rightEndBitLine         = plot(     outputStruct.rightEndBit.xBeg, ...
                                  	outputStruct.rightEndBit.yBeg, 'r.-');

% plot the points that a subject to failure (local maxima)                                
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
xlabel('x (m)')
ylabel('y (m)')

if longPlotOn

% b) displacmement (shear, normal)
subplot(numPlots, 1, subPlotCount) 
subPlotCount = subPlotCount+1;
hold on 

yyaxis left
shearDisplacementPlot   = plot(outputStruct.faultInfo.xBeg, ...
                               outputStruct.faultStress.DS);
ylabel('Shear Displacement (m)')
                           
yyaxis right
normalDisplacementPlot  = plot(outputStruct.faultInfo.xBeg, ...
                               outputStruct.faultStress.DN);
ylabel('Normal Displacement (m)')
xlabel('x (m)')
                           
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
xlabel('x (m)')
ylabel('(MPa)')

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
xlabel('x (m)')
ylabel('(MPa)')

end
    % nested function: 
   	function plotGROWOutput(what2plot)
            % plottting grow output
            % input:
            % what2plot: cell array, or string with thing to plot
            %   'geometry': plot crack geometry
            
            
            if strcmp(inputParameters.user.runGROWModel, 'yes') % no error meassage here because it should be issued in earlier functions
                if ~strcmp(inputParameters.user.plotGROWModel, 'yes')
                    if strcmp(inputParameters.user.plotGROWModel, 'no')
                        disp('Grow output not plotted, change inputParameters.user.plotGROWModel to: ''yes''')
                    else
                        error('''plotGROWMOdel'' must be ''yes'' or ''no''')
                    end
                    return
                end
                
                if any(strcmp(what2plot,'geometry'))
                    hold on
                    for iSegment = 1:(size(crackLocations,1)+1)
                        
                        xBeg = GROWOutputStruct.segments(iSegment).xbeg;
                        yBeg = GROWOutputStruct.segments(iSegment).ybeg;
                        xEnd = GROWOutputStruct.segments(iSegment).xend;
                        yEnd = GROWOutputStruct.segments(iSegment).yend;
                        
                        dx = xEnd - xBeg;
                        dy = yEnd - yBeg;
                        
                        segmentLine = quiver(xBeg, yBeg, dx, dy);
                        segmentLine.ShowArrowHead   = 'off';
                        segmentLine.AutoScale       = 'off';
                        segmentLine.Marker          = '.';
                      
                    end % for segments
                else
                    error('input to function plotGROWOutput must be one, or many, of ''geometry''')
                end % if geomtry is called
            end % if grow model is run
    end % function to plot grow output
    
end % function to plot output

% 8) organize files to reduce clutter in directory
    function organizefiles
    
% make new directory for the files from this run

% make an output directory if it does not already exist
outputDirectoryName = 'output';
if ~isdir(outputDirectoryName)
    system(['mkdir ', outputDirectoryName])
end
    
newDirectoryName    = [outputDirectoryName, '/', fileName];

if ~isdir(newDirectoryName)
command             = sprintf('mkdir %s', newDirectoryName);
system(command);
end 

% move all file that start with 'fileName' into new directory
command = sprintf('mv %s %s', [fileName,'.*'], [newDirectoryName,'/']);
system(command);
command = sprintf('mv %s %s', 'grow_specs*', [newDirectoryName,'/']);
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
 
% save matlab content including output figure and input specs
savefig([newDirectoryName,'/',fileName]);
save([newDirectoryName,'/',fileName], 'inputParameters', 'GROWInputParameters');

system('rm end_bits.txt');
system('rm input');
system('rm spliced_profile_temp.txt');
system('rm fault_header.txt')
profileName = [fileName,'_profile.txt']; system(['rm ',profileName]);

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------- User Input -----------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [inputParameters, GROWInputParameters] = userInput()
    inputParameters = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % THIS IS LIKELY THE ONLY SECTION A USER MAY WANT TO CHANGE %%%%%%%%%%%
    % material parameters:
    
    E                       = 3000;     % Young's Modulus MPa
    nu                      = 0.25;     % Poisson's ratio
    T                       = 10;       % Tensile strength of the host rock (MPa) | 5-10
    S_0                     = 25;        % Shear strength of the host rock (MPa) | 25-60 (realistic amounts)
    static_friction         = 0.6;      % static friction of crack elements
    dynamic_friction        = 0;        % dynamic friction of crack elements
    critical_slip_distance  = 10^-5;    % critical slip distance (m)
    shear_stiffness         = 10^10;    % shear stiffness 
    normal_stiffness        = 10^10;    % normal stiffness
    sliding_cohesion        = 0;        % sliding cohesion (resistance of fault elements to tension)
    
    % boundary conditions
    
    lithostatic_stress      = -40; 
    shear_stress            = -40;
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MORE TECHNICAL CHANGES:
    
    % user choices

    inputParameters.user.showFric2DOutput          = 'yes';     % plot output from fric2d ('yes' or 'no')
    inputParameters.user.cracks2grow               = 'max coulomb';
        % options:
        % 'all':           all the coordinates specified in input (tensile+shear)
        % 'tensile':       only cracks that failed in tensile stress (listed in
        %                  tensileFailureCoord)
        % 'shear'          only cracks that failed in shear stress (listed in
        %                  shearFailureCoord)
        % 'max tensile'    most tensile crack
        % 'max coulomb'    most sheared crack  
        % 'max tensile and max coulomb' both most shear and most tensile
                                                                 
    inputParameters.user.runGROWModel              = 'yes' ;     % run grow ('yes' or 'no')
    inputParameters.user.plotGROWModel             = 'yes' ;     % plot grow output 'yes' or 'no' (make sure runGROWmodel is set to yes as well)
    
    % GROW inputs
    
    GROWInputParameters = [];
    
    GROWInputParameters.growFileName = 'grow_specs';
    
    % grow input  parameters parameter
   
    % this will be part of the command to run GROW (Must be digits - no
    % decimal points)
    GROWInputParameters.angleResolution         = 15;
    GROWInputParameters.startAngle              = 90;
    GROWInputParameters.endAngle                = 270;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STORING INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    inputParameters.material.E                     = E;                     % Young's Modulus MPa
    inputParameters.material.nu                    = nu;                    % Poisson's ratio

    inputParameters.material.c                     = T;                     % cohesion (MPa)
    inputParameters.material.phi                   = atand(static_friction); % angle of intenral friction  (degrees)

    % boundary conditions

    inputParameters.boundaryConditions.sigxx       = lithostatic_stress;   	% normal stress parallel to faults MPa
    inputParameters.boundaryConditions.sigyy       = lithostatic_stress;    % normal stress perpendicular to faults MPa
    inputParameters.boundaryConditions.shear       = shear_stress;          % shear stress on faults  MPa

    % parameters tacked on to input file
    fault_flaw              = [];
    fault_flaw.stiffS       = shear_stiffness;          % shear stiffness
    fault_flaw.stiffN       = normal_stiffness;         % normal stiffness
    fault_flaw.T            = T;                        % Tensional strenght
    fault_flaw.S_0          = S_0;                      % Internal shear strength
    fault_flaw.cohes        = sliding_cohesion;       	% sliding cohesion
    fault_flaw.friction_s   = static_friction;          % static friction
    fault_flaw.friction_d   = dynamic_friction;       	% dynamic friction
    fault_flaw.L            = critical_slip_distance; 	% critical sliding distance
    
    GROWInputParameters.fault_flaw = fault_flaw;
    
    % software:
    GROWInputParameters.GROW_perl_fileName = 'GROW_OCT_2017.pl';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------- section functions ----------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input
function [X,Y, plotOption] = checkinput(fileName, X,Y, input, numInput)
        
        if isa(fileName, 'string')
            error('input '' fileName must be of ype ''string''');
        end
    
        if numInput < 3
            error('There must be at least three inputs: an X array, a Y array and a string specifying the file name')
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
         
         if ~isempty(input)
            plotOption = input{1};
         else
             plotOption = 'long';
         end
         
end

% check whether necessary files are in folder:
function check_directory(GROWInputParameters)

% check directory

checkFile(GROWInputParameters.GROW_perl_fileName);
checkFile('simple_shear_bc.pl');
checkFile('Wext.pl');

% compile fric2d if hasnt been done so already and move file into working
% directory
if (exist('fric2d', 'file') == 0)
    
    % IMPORTANT: if fric
    path2fric2dCompiler = 'GROW/fric2d_source_code';
    
    disp('fric2d is not in the working directory, we will try to take care of it...')

    currentFolder = pwd;
    
    try %to compile fric 2d
    cd(path2fric2dCompiler);
    system('make');
    cd(currentFolder);
    system(['cp ', path2fric2dCompiler, '/fric2d fric2d']);
    disp('got it')
    catch
        error(['the fric2d executable is not in working directory, a compiler was not found in the following directory: ', ...
               path2fric2dCompiler]);
    end
end

    function checkFile(fileName)
        
        if exist(fileName, 'file') ~= 2
            errorMessage = [fileName, ' is not in directory, you need it to be'];
            error(errorMessage)
        end  
    end


end

%% 1) make the input file for the fric2D program
function [meanSegmentLength, X] = makeinputfile(fileName,X,Y, inputParameters)

% create profile: 
profileName                 = [fileName,'_profile.txt']; % dont change this same name used in 6b

padding = 2;
boxAspectRatio = 0.3; % height to length ratio of the boundary box

[X, meanSegmentLength,boxSize] = makeprofile(profileName,X,Y,padding, boxAspectRatio);

% create boundary condition (creates a file named 'input')
    % input ptspacing based on number of points
    % simple_shear_bc.pl 0.06 0.15 0.0002
    % (based on profile1)

    input_section_fileName = 'input';
    
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

% calculate number of endbit segments on either side of input fault
% geometry (magic number 4 is used to calulate half the distance from the
% tip of the fault geometry to the end  side of the box)

numBufferSegments = floor(boxSize(1)/(padding*4*meanSegmentLength));
make_end_bits(X,Y, numBufferSegments, meanSegmentLength, boxSize, inputParameters);

% add tag line for the first fault
fault_header_fileName = 'fault_header.txt';
fid=fopen(fault_header_fileName,'w');
fprintf(fid, 'fault segment_0 yes no no\n');
fclose(fid);

% concatenate to the input file
command = sprintf('cat %s %s %s %s > %s',   ...
                  input_section_fileName,   ...
                  fault_header_fileName,    ...
                  profileName,              ...
                  'end_bits.txt',           ...
                  [fileName, '.in']         );
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
MEAN_SEGMENT_LENGTH = round(MEAN_SEGMENT_LENGTH,2,'significant');

% boxSize = padding.*[RX,RY]; % not great because the box size changes as a
% function of topography, which is bad if set in a loop

% instead:

boxSize = zeros(1,2); 
boxSize(1) = padding(1)*RX;
boxSize(2) = boxSize(1)*boxAspectRatio;

% round boxsize to nearest multiple of boundary element length
boxSize = ceil(boxSize/MEAN_SEGMENT_LENGTH)*MEAN_SEGMENT_LENGTH;

% errors/warnings related to profile in relation to box size
if max(abs(Y))>boxSize(2)*0.1
    if max(abs(Y))>boxSize(2)/2
        error('Profile topography is too large and intersect the boundary box, choose a larger aspect ratio for the box by changing the vatiable ''boxAspectRatio'' ')
    end
    warning('Profile topography is large, consider choosing a larger aspect ratio for the box by changing the vatiable ''boxAspectRatio'' ');
end

% round to increments of the length of elements

newX = X - min(X) + boxSize(1)/2 - RX/2; % shift X to be centered in the observation box

fid=fopen(fileName,'w');

for n = 1:(length(newX)-1)
    fprintf(fid, ...
        '1 %f %f %f %f 0 1.00E+10 0 0 0 0.6 0 0\n', ...
        [newX(n),  Y(n),  newX(n+1),   Y(n+1)]');
end

fprintf(fid, '\n');
fclose(fid);

end

% 1b) add end bits to input file according to numBuffer segments (figiting
% needed)

    function [] = make_end_bits(X,Y,numBufferSegments, meanSegmentLength, BOXSIZE, inputParameters) % MAKE ThE FUNCTION CALL MATCh

    intro =    {'*FRACTURES'; ...
    '*num	xbeg	ybeg	xend	yend	kode	bvs     bvn     bvs     bvn'; ...
    '*----	----	----	----	----	----	------	----	-----	--------'};
    
    tempFileName = 'end_bits.txt';
    fileID = fopen(tempFileName,'w');
    fprintf(fileID,'%s\r\n',intro{:});
    fclose(fileID);
    % find ends of the fault
    X1      = X(1);
    Y1      = Y(1);
    Xend    = X(end);
    Yend    = Y(end);

    % find end of each end bit
    start1  = [X1-numBufferSegments*meanSegmentLength,Y1];
    end1    = [X1,Y1];
    start2  = [Xend,Yend];
    end2    = [Xend+numBufferSegments*meanSegmentLength,Yend];
    
    modelHeight = BOXSIZE(2); % CHECK THIS
    make_one_EB('left' ,tempFileName, start1,end1,numBufferSegments,meanSegmentLength, modelHeight, inputParameters);
    make_one_EB('right',tempFileName, start2,end2,numBufferSegments,meanSegmentLength, modelHeight, inputParameters);

        function [] = make_one_EB(EBName, fileName, EBStart, EBEnd,numBufferSegments, meanSegmentLength, MODELHEIGHT, inputParameters)
            
            % function to make one end bit to avoid redundancy...
            
            % make set of X coordinates for each end bit
            % Left end bit (LEB):
            EBX      = linspace(EBStart(1),EBEnd(1),numBufferSegments+1);
            
            % make start and end
            EBXBeg  = EBX(1:end-1);
            EBXEnd  = EBX(2:end);
                        
            EBY      = linspace(EBStart(2),EBEnd(2),numBufferSegments+1);
            
            % make start and end
            EBYBeg  = EBY(1:end-1);
            EBYEnd  = EBY(2:end);
            
            % ... and midpoints
            EBXMid  = (EBXBeg+EBXEnd)/2;
            EBYMid  = (EBYBeg+EBYEnd)/2;
                                 
            EBUs    = calc_end_bit_displacement(EBXMid,EBYMid, ...
                                                MODELHEIGHT, ...
                                                inputParameters.material.E  , ...
                                                inputParameters.material.nu , ...
                                                inputParameters.boundaryConditions.sigxx, ...
                                                inputParameters.boundaryConditions.sigyy, ...
                                                inputParameters.boundaryConditions.shear, ...
                                                meanSegmentLength);
            
            % header section
            FID = fopen(fileName,'a');
            fprintf(FID, 'fracture %s no no\n',EBName);
            kode = 3;
            bvs = 0;
            bvn = 0;
            
            % individual elements
            for iElement = 1:numBufferSegments               
                
                ROW1    = [ 1, ...
                            EBXBeg(iElement), ...
                            EBYBeg(iElement), ...
                            EBXEnd(iElement), ...
                            EBYEnd(iElement), ...
                            kode, ...
                            bvs,...
                            bvn,...
                            EBUs(iElement),...
                            inputParameters.boundaryConditions.sigxx ];
                        
                fprintf(FID, '%f %f %f %f %f %f %f %f %f %f\n',ROW1);
                
            end
            fprintf(FID, '\n');
            fclose(FID);
            
        end

        function [US] = calc_end_bit_displacement(x,y,h,E,nu,sigyy,sigxx,shear,elementSize)
            % analytical solution for the displacement (u_s) on the end bits  

            % Input (1 by 1 or 1 by N)
            % [x,y]:            center coordinate of the endbit element
            % E:                Yougn's Mofulous (MPa) 1 by 
            % nu:               Poisson's ratio
            % [sigxx, sigyy]:   Lithostatic stress (MPa)
            % shear:            Shear stress (MPa)
            % elementSize:      Size of boundary elements (2x slip) (m)

            % analytical calculation
            US = -(x./E) .* ((1-nu.^2).*-sigxx - nu.*(1+nu).*-sigyy) +...
                 ( -shear  .* (y+h/2)   .*2*(1+nu)./E) - ...
                 elementSize/4;
             
        end % function [US] = calc_end_bit_displacement(x,y,E,nu,sigyy,sigxx,shear,elementSize)
    end %  function [] = make_end_bits(X,Y,numBufferSegments, BVSTop, BVSBot)

end

%% 2) run fric2D
function runfric2D(fileName)

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

dataBlocks = read_blocks_of_numerical_data(outputFileName, 10);
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
                               'kode'           , ...
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
    
    % DOING GEOLOGY FOR A MOMENT HERE
    sigma1      = -eigenValues(2);
    sigma3      = -eigenValues(1);
    
    tau_m       = (sigma1-sigma3)/2;
    sigma_m     = (sigma1+sigma3)/2;
    
    coulombCriterion(iBE) = -tau_m + (sigma_m*sind(phi)+c*cosd(phi));
    tensileCriterion(iBE) = -sigma3-c;
 
end

% find points of failure (findpeaks is used to find local maxima in case
% there are segments of points that are above failure)
tempCC = coulombCriterion;
tempTC = tensileCriterion;

%tempCC(tempCC<0) = 0; % remove section that are not going to fail
%tempTC(tempTC<0) = 0;

[~,shearFailureInd  ] = findpeaks(tempCC,'SortStr', 'descend'); % find peak stresses and sort their indices
[~,tensileFailureInd] = findpeaks(tempTC,'SortStr', 'descend');

shearFailureCoord     =  [X(shearFailureInd+1), Y(shearFailureInd+1)];                           
tensileFailureCoord   =  [X(tensileFailureInd+1), Y(tensileFailureInd+1)];                     
                     
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

CRACKLOCATION = [];
if      strcmp(userChoice, 'all');          try CRACKLOCATION = [shearFailureCoord;tensileFailureCoord];    end 
elseif  strcmp(userChoice, 'tensile');      try CRACKLOCATION = tensileFailureCoord;                        end
elseif  strcmp(userChoice, 'max tensile');  try CRACKLOCATION = tensileFailureCoord(1,:);                   end % max is the first location because the list is sorted   
elseif  strcmp(userChoice, 'shear');        try CRACKLOCATION = shearFailureCoord;                          end  
elseif  strcmp(userChoice, 'max coulomb');  try CRACKLOCATION = shearFailureCoord(1,:);                     end % max is the first location because the list is sorted
elseif  strcmp(userChoice, 'max tensile and max coulomb') 
    try CRACKLOCATION = [tensileFailureCoord(1,:);shearFailureCoord(1,:)];                                  end 
else;                                       error('choosecracks ''userChoice'' input must be one of: ''all'', ''tensile'',''shear'',''max tensile'' or ''max coulomb''')
end

% note that this is effectively the 'catch' to the try's above
if isempty(CRACKLOCATION) 
    disp('No cracks will fail given the specified conditions')
    return
end

end

%% 5b) use outputs from fric2D and failure asseessment to track fracture growth
function [] = runGROW(X, Y, crackLocations, inputParameters, GROWInputParameters)
% this function runs GROW (GROwth by Work-minimization)

% input:
% fileName: string with the name of a fric2D input file 

GROWInputFileName = [GROWInputParameters.growFileName, '.in'];

fault_flaw  = GROWInputParameters.fault_flaw;

intact_flaw.xcoord      = crackLocations(:,1);
intact_flaw.ycoord      = crackLocations(:,2);

% check if the function was called
if ~strcmp(inputParameters.user.runGROWModel, 'yes')
    if strcmp(inputParameters.user.runGROWModel, 'no')
        disp('No GROW Model was lauched, change ''inputParameters.user.runGROWModel'' input to ''yes'' if so desired')
    else 
        error('inputParameters.user.runGROWModel in input section must be ''yes'' or ''no''')
    end
    return

end

% split fault into two and add *Crack tag to the on of the fauts. One of
% the faults will be allowed to grow from one of its ends

% run for-loop over fault block
% when hit a failure point add header line
% \n\n
% fault rough-n yes yes no
% *Tag stiffS stiffN T S_0 cohes friction-s friction-d L 
% *Crack  1e10 1e10 1.000000 40.000000 0.000000 0.600000 0.000000 0.000010

% the above are input parameters (make these editable)

% 

% cat together input section, boundary conditions, spliced faults, endbit


done = 1;
if done
    
    numCrack = size(crackLocations,1);
    
    % splice profile:
    % this is not exactly the smartest way to do this...
    spliced_profile_fileName = 'spliced_profile_temp.txt';
    fid = fopen(spliced_profile_fileName,'w');
    
    % print first header line
    fprintf(fid, 'fault segment_0 yes no yes\n');
    fprintf(fid, '*%s %s %s %s %s %s %s %s %s\n', ...
                             'Tag'                        , ...
                             'stiffS'                       , ...
                             'stiffN'                       , ...
                             'T'                            , ...
                             'S_0'                          , ...
                             'cohes'                        , ...
                             'friction-s'                   , ...
                             'friction-d'                   , ...
                             'L'                            );
                         
                fprintf(fid, '*%s %f %f %f %f %f %f %f %f\n', ...
                                'Crack'                     , ...
                                fault_flaw.stiffS           , ...
                                fault_flaw.stiffN           , ...
                                fault_flaw.T                , ...
                                fault_flaw.S_0              , ...
                                fault_flaw.cohes            , ...
                                fault_flaw.friction_s       , ...
                                fault_flaw.friction_d       , ...
                                fault_flaw.L                );
    
    crackCount = 0;
    for n = 1:(length(X)-1)
        
        if any(X(n) == crackLocations(:,1))
            
            crackCount = crackCount+1;
            
            if crackCount ~= numCrack

                % printing the following to the file:         

                %             fault   rough1    yes     no      yes    
                % *Tag stiffS stiffN T S_0 cohes friction-s friction-d L 
                % *Crack  1e10 1e10 1.000000 1.000000 0.000000 0.600000 0.000000 0.000010

                fprintf(fid, '\n');
                fprintf(fid, 'fault segment_%i yes no yes\n'  , crackCount);
                fprintf(fid, '*%s %s %s %s %s %s %s %s %s\n', ...
                             'Tag'                        , ...
                             'stiffS'                       , ...
                             'stiffN'                       , ...
                             'T'                            , ...
                             'S_0'                          , ...
                             'cohes'                        , ...
                             'friction-s'                   , ...
                             'friction-d'                   , ...
                             'L'                            );
                         
                fprintf(fid, '*%s %f %f %f %f %f %f %f %f\n', ...
                                'Crack'                     , ...
                                fault_flaw.stiffS           , ...
                                fault_flaw.stiffN           , ...
                                fault_flaw.T                , ...
                                fault_flaw.S_0              , ...
                                fault_flaw.cohes            , ...
                                fault_flaw.friction_s       , ...
                                fault_flaw.friction_d       , ...
                                fault_flaw.L                );
                             
            else % if crackCount ~= numCrack
                
                % if it is the last crack, we do not want to grow the end
                % of the fault
                
                fprintf(fid,'\n');
                fprintf(fid,'fault segment_%i yes no no\n', crackCount);
            end
        end
        
        fprintf(fid, ...
            '1 %f %f %f %f 0 1.00E+10 0 0 0 0.6 0 0\n', ...
            [X(n),  Y(n),  X(n+1),   Y(n+1)]');
        

    end
    
    fprintf(fid,'\n');
    fclose(fid);
    
    % concatenate the input file (generate by simple_shear.pl), the pro
    command = sprintf(  'cat %s %s %s > %s', ...
                        'input'                             , ...
                        spliced_profile_fileName            , ...
                        'end_bits.txt'                      , ...
                        GROWInputFileName);
    system(command);
    
    
    command = sprintf('perl %s %s %.0f %.0f %.0f'              , ...
        GROWInputParameters.GROW_perl_fileName, ...
        GROWInputFileName                   , ...
        GROWInputParameters.angleResolution , ...
        GROWInputParameters.startAngle    	, ...
        GROWInputParameters.endAngle    	);
    
    system(command);
    
    
end

% stuff from Michele
%- add lines in the Fault conditions section: (make sure parameters
% are rigth for sandstone)

% *Tag fault xcoor ycoor length grow_both? stiffS stiffN T S_0 cohes friction-s friction-d L
% *Flaw-Fault 			 				   	   0	1e10		0	0	0	1	     	1			0.0	
% *Flaw-Intact 	rough	0.055424 0.01173379	0.0006	yes         1e10	  1e10		5	50	0	0.6	     	0			1e-5

% to run: ./GROW_JULY_2016.pl bump.in 45 45 315


end

%% 6) extract Grow output 
function GROWOutputStruct = extractGROWoutput(crackLoc, growFileName, inputParameters)

% extracts numeric blocks of information out of GROW output file and saves
% them into matlab structure with fields names corresponding (sometimes
% with small adjustments to header lines)

% % input: % fileName: name of the file in a string without the extension
% (here % assumed to be  .out). IMPORTANT: output file is assumed to have 9
% data % blocks + x blocks per flaw. Namely 2 boundary line blocks,
% followed by 4 'end bit' blocks and % 3 fault info blocks. The current
% state does not extract 'end bits' % information.

% output:
% GROWOutputStruct: structure with with fields corresponding to blocks of data
%   segments: 1 by N (N is numuber of fault segments, or, equivalently, the
%   number of defects+1) structure with fields corresponding to .eff
%   headers for the crack information (num    xbeg    ybeg    xend    yend
%   stiffS stiffN ten_str init_coh slid_coh stat_fric dy_fric crit_slip)
% 
% POSSIBLY MORE COMING...

% function calls:
% read_blocks_of_numeric_data
% assign2struct

% usage:
% extractoutput('sample_section_test')

% 2017-04-17: tested and functions properly using the file
% sample_section_test.out


% check if this function should be run:
if ~strcmp(inputParameters.user.runGROWModel, 'yes')
    if strcmp(inputParameters.user.runGROWModel, 'no')
        disp('No GROW Model was lauched, change ''inputParameters.user.runGROWModel'' input to ''yes'' if so desired')
    else 
        error('inputParameters.user.runGROWModel in input section must be ''yes'' or ''no''')
    end
    return
end
        
outputFileName = [growFileName,'.eff'];

% load in block of data (TERRIBLE IDEA BUT FUCK IT): 

% determine howmany cracks where grown
crackLocDim = size(crackLoc);
numCracks = crackLocDim(1);
crackBlocks = numCracks+1; % make sure this is correct
boundaryBlocks = 1;
endBitBlocks = 2;
expectedFric2dDataBlocks = boundaryBlocks+endBitBlocks; % this is here to make sure that I'm not fucking up the parsing of the output


expectedDataBlocks = expectedFric2dDataBlocks + crackBlocks;

% extracts blocks of numeric data
dataBlocks = read_blocks_of_numerical_data(outputFileName, 10,' |\t'); % 10 is the minimum length of string
numBlocks  = length(dataBlocks);

% double checking that the number of blocks is right
if numBlocks ~= expectedDataBlocks
    error('the GROW output is not of the expected size consider double-checking output parsing')
end

% Get crack info:
numSegments = numCracks+1;

% initialize structure
segmentHeaders = {  'num'               , ...
                    'xbeg'              , ...
                    'ybeg'              , ...
                    'xend'              , ...
                    'yend'              , ...
                    'stiffS'            , ...
                    'stiffN'            , ...
                    'tens_str'          , ...
                    'init_coh'          , ...
                    'stat_fric'         , ...
                    'dy_fric'           , ...
                    'crit_slip'         };
blankData = zeros(1,length(segmentHeaders));
segmentStruct = assign2struct(blankData,segmentHeaders);
segmentStruct(numSegments).(segmentHeaders{1}) = 0;         

for iSegment = 1:numSegments
    % exctract data for each crack
    
    % steps:
    % define headers (what will they be?)
    
    % get block of data
    % get the geometry
    
    % What I need to do here:
    % figure out how many blocks of data will be produced from GROW
    
    % for now the key blocks of data that I want are the crack geometry
    % need to know what the header sections will be as well
    
    segmentBlock = dataBlocks(boundaryBlocks+iSegment);
    segmentStruct(iSegment) = assign2struct(segmentBlock{:},segmentHeaders);

end % for iCrack = 1:numCrack

% output is stroredin mater structure 
% not that for now this is a bit odd,m but I am leave the possibility for
% extracting other fields from the efficient file (or others, for that
% matter)
GROWOutputStruct = [];
GROWOutputStruct.segments = segmentStruct; 
    
end 

%% 7) visually diaplay output information (function is located in embeded function section)

%% 8) Organize files (function is located in embeded function section)

%% 9) check model - issue potential warnings

function checkmodel(X,Y,inputParameters, outputStruct, GROWOutputStruct)

% displacement limitations
elementLength = sqrt((X(2)-X(1))^2 + (Y(2)-Y(1))^2);
messageOut = sprintf('Element length: %f, Max displacement: %f', ...
                     elementLength, ...
                     max(abs(outputStruct.faultStress.DS)));
disp(messageOut)

   if any(abs(outputStruct.faultStress.DS) > 0.5*elementLength)
        warning('displacement on individual elements is too large, consider making elements longer or reducing stresses')
   end
   
end

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


