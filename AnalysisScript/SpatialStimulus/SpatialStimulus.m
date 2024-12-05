% Analysis for spatial stimuli

% Written by Felice Dunn

%% Getting the full responses and frame rates

% change the number of the first children to move from one cell to the next
epoch = tree.children(1).children(1).children(1).children(1).children(1).epochList 

data = riekesuite.getResponseMatrix(epoch, 'Amplifier_Ch1');

%% Save stimuli in the same folder 

clear allSeeds 
  
[a, b] = size(data);
   
for n = 1:a
    allSeeds(n) = epoch.elements(n).protocolSettings.get('noiseSeed')
end
 
format long, allSeeds

 
lines = 83; % 20200706 new conversion factor gives 83 lines of 40um bars
frames = 52*30;
mn = 0.5; % mean
contrast = 0.3; % sd
binary = 1; % binary or not


saveStims(lines,frames,mn,contrast,binary,allSeeds)
