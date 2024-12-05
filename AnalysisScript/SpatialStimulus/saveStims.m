% Written by Felice Dunn

function saveStims(lines,frames,mn,contrast,binary,allSeeds)
numStims=length(allSeeds);
for i=1:numStims
    thisSeed=allSeeds(i);
    obj.noiseStream = RandStream('mt19937ar', 'Seed',thisSeed);
    if binary
        stim=(255.*2*mn.* (obj.noiseStream.rand(lines,frames)>0.5));
    else
        stim=(255.*(mn + (mn*contrast) * (obj.noiseStream.randn(lines,frames))));
    end
    
    suffix=i-1;
    stimStr = ['stim', num2str(suffix),'.txt'];
    
    save(stimStr','stim','-ASCII')
end





%%% as of November 2016, the command should look like this for 50 seconds of 40um lines:

% saveStims(35,1505,0.5,0.3,1,allSeeds_1201c1)

% with the appropriate allSeeds_(date)(cell#) containing the noise seeds in order as presented in the matlab tree gui


