function sample_run(varargin)

% name of the run 

generalRunName = 'sample_run';

% make the desired profile

numPts = 40;
profileType =  'single_asperity';

numberOfRuns = 1;

if strcmp(profileType,'single_asperity')
    
    % in construction
    
    asperityType            = {'sine'};     % type of asperity (sine, traingle, box, step, etc - see the avialble functionaltity below)
    absoluteAsperityLength  = 0.01;         % width of the asperity in m
    absoluteAsperityHeight  = 0.004;        % height of the asperity in m
    absoluteFaultLength     = 0.07;          % fault length in m
    %padding                 = 3;           % DEPRECATED : seems like a bad idea to change the box size (function of fault length) % padding on either side of the asperity that will just be flat (a factor of the asperityLength

    % duplicate input (make another specific array afterward if so desired
    asperityTypeArray            = repmat(asperityType,              1,numberOfRuns); % type of asperity (sine, traingle, box, step, etc - see the avialble functionaltity below)
    absoluteAsperityLenghtArray  = repmat(absoluteAsperityLength,    1,numberOfRuns); % width of the asperity in m
    absoluteAsperityHeightArray  = repmat(absoluteAsperityHeight,    1,numberOfRuns); % height of the asperity in m
    %paddingArray                 = repmat(padding,                   1,numberOfRuns); % padding on either side of the asperity that will just be flat (a factor of the asperityLength
   
    absoluteFaultLengthArray     = repmat(absoluteFaultLength,       1,numberOfRuns); 
    
end

numPtsArray = repmat(numPts, 1,numberOfRuns);

%%%%% make custom array of length Nunmber of runs here: %%%%%%%%%%%%%%%%%%%
% e.g.:
%numPtsArray = ceil(linspace(10,300,numberOfRuns));
% absoluteAsperityLenghtArray  = 0.01./(1:numberOfRuns);
% absoluteAsperityHeightArray = 0.001./(1:numberOfRuns);
% asperityTypeArray            = {'sine','triangle', 'box','step','jog'}; if length(asperityType ~= numberOfRuns); error(asperityTypes must be of length number of runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iRun = 1:numberOfRuns
    
if strcmp(profileType,'single_asperity')
    
    asperityType            = string(asperityTypeArray{iRun});
    absoluteAsperityLength  = absoluteAsperityLenghtArray(iRun);
    absoluteAsperityHeight  = absoluteAsperityHeightArray(iRun);
    % padding                 = paddingArray(iRun);     
    numPts                  = numPtsArray(iRun);
    absoluteFaultLength     = absoluteFaultLengthArray(iRun);
    
 
    
    pointSpacing = absoluteFaultLength/numPts; % spacing of points in 
   
    X = (0:numPts-1)*pointSpacing;
    
    numAsperityPts = ceil(numPts/absoluteFaultLength*absoluteAsperityLength);
    startInd = floor(numPts/2-numAsperityPts/2);
    
    if strcmp(asperityType, 'sine') % Hann window
        Y = [zeros(1,startInd), ...
                tukeywin(numAsperityPts,1)'*absoluteAsperityHeight, ...
                zeros(1,numPts-startInd-numAsperityPts)];
    
    elseif strcmp(asperityType, 'triangle') % bartlet window
        Y = [zeros(1,startInd), ...
            window(@bartlett, numAsperityPts)'*absoluteAsperityHeight, ...
            zeros(1,numPts-startInd-numAsperityPts)];
        
    elseif strcmp(asperityType, 'box') % tuckyWin
       Y = [zeros(1,startInd), ...
                tukeywin(numAsperityPts,0)'*absoluteAsperityHeight, ...
                zeros(1,numPts-startInd-numAsperityPts)];
    
    elseif strcmp(asperityType, 'step') % step
       Y = [zeros(1,floor(numPts/2)), ...
              ones(1,ceil(numPts/2))*absoluteAsperityHeight];
    elseif strcmp(asperityType, 'jog') % diagonal step
       Y = [zeros(1,startInd), absoluteAsperityHeight/numAsperityPts*(1:numAsperityPts), ...
           ones(1,numPts-startInd-numAsperityPts)*absoluteAsperityHeight];
    end
    
elseif strcmp(profileType,'real_profile')
    % to do
end
    
fileName = sprintf('%s_%s_%s_%s_run%i',generalRunName, profileType, asperityType, date, iRun);
                
% run the work flow
tic
fric2d_workflow(fileName,X,Y)
runTime = toc;
sprintf('run time was: %f seconds', runTime)
end


end
