% name of the run 

generalRunName = 'sample_run';

% make the desired profile

numPts = 100;
profileType =  'single_asperity';

if strcmp(profileType,'single_asperity')
    
    asperityType = 'sine'; % type of asperity (sine, traingle, box, step, etc - see the avialble functionaltity below)
    absoluteAsperityLenght = 0.01; % width of the asperity in m
    absoluteAsperityHeight = 0.001; % height of the asperity in m
    padding   = 3; % padding on either side of the asperity that will just be flat (a factor of the asperityLength

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pointSpacing = (2*padding+1)*absoluteAsperityLenght/numPts; % spacing of points in 
   
    X = (0:numPts-1)*pointSpacing;
    
    numAsperityPts = ceil(numPts/(2*padding+1));
    startInd = floor(numPts/2-numAsperityPts/2);
    
    if strcmp(asperityType, 'sine') % Hann window
        Y = [zeros(1,startInd), ...
                tukeywin(numAsperityPts,1)'*absoluteAsperityHeight, ...
                zeros(1,numPts-startInd-numAsperityPts)];
    
    elseif strcmp(asperityType, 'triangle') % bartlet window
        Y = [zeros(1,startInd), ...
            window(@bartlet, numAsperityPts)'*absoluteAsperityHeight, ...
            zeros(1,numPts-startInd-numAsperityPts)];
        
    elseif strcmp(asperityType, 'box') % tuckyWin
       Y = [zeros(1,startInd), ...
                tukeywin(numAsperityPts,1)'*absoluteAsperityHeight, ...
                zeros(1,numPts-startInd-numAsperityPts)];
    
    elseif strcmp(asperityType, 'step') % step
       Y = [zeros(1,floor(numPts/2)), ...
              ones(1,ceil(numPts/2))*absoluteAsperityHeight];
    end
    
elseif strcmp(profileType,'real_profile')
    % to do
end
    
fileName = [generalRunName,'_', ...
                    profileType, '_', ...
                    asperityType, '_', ...
                    date];
                
% run the work flow
fric2D(fileName,X,Y)
