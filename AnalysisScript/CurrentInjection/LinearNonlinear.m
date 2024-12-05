% Linear-Nonlinear Model for Spikes

%% To calculate the spike number at each node separately
 
epoch = tree.children(1).children(1).children(1).children(1).epochList;
 
data = riekesuite.getResponseMatrix(epoch, 'Amplifier_Ch1');
 
spikes = SpikeDetector_SK(data);

% To get raster plot for spikes (run only if needed for the figure)

% From 'spike_analysis_script_JL/RC.m'

% To get raster plots, run above three lines for the node you want spikes
% for, then navigate to the folder you want to save them in and run:
 
spiketrace = zeros(1,18000);  % change based on data size
% spiketrace(spikes.sp{1}) = 1; % number in curly brackets is which trace of the 4 or 5 at that flash strength
spiketrace(spikes.sp) = 1;  % don't need {} when there's only 1 epoch
spiketrace = spiketrace';

figure
plot (spiketrace);

%% To get the stimulus vector for injection

condition = 'Amplifier_Ch1'; % change the condition to get the light stimulus

stimulus = riekesuite.getStimulusMatrix(epoch, condition); 

% 20220714: use getStimulusMatrix and getStimlusVector from +riekesuite
% folder not Juan's. 

%% spike analysis faster

num = length(spikes.sp);
STA = [];
spiketriggered = [];
STA_mean = [];

clear STA_mean

timebefore = 1000; % number of points to draw before the spike (5000) 1000 for current injection
for n = 1:num
    clear points cutoff point_cutoff spiketriggered
    points = spikes.sp{n};
    cutoff = find(points> timebefore);% remove = because ran into a problem when = timebefore
    point_cutoff = points(cutoff)+8; %throw away stimulus that is too short from the timebefore (added 20 to account for delay in the stimulus 20180125 figure out what this should be)
    ave = mean(stimulus(n,:));
    spiketriggered = zeros(1, timebefore + 1);
    if isempty(cutoff) == 0
        for cnt = 1:length(point_cutoff);
            spiketriggered = spiketriggered + ((stimulus(n, (point_cutoff(cnt)-timebefore):point_cutoff(cnt)))); % subtract off mean of stimulus to zero
        end
        STA(n, :) = (spiketriggered)./cnt;
    end
end

STA_mean = mean(STA);

figure
plot(STA_mean, 'k')

%% Sccount for late spikes

num = length(spikes.sp);
STA = [];
spiketriggered = [];
STA_mean = [];
buffer = 20;

clear STA_mean

timebefore = 1000; % number of points to draw before the spike (5000) 1000 for current injection
for n = 1:num
    clear points cutoff point_cutoff spiketriggered
    last = length(stimulus(n,:));
    points = spikes.sp{n};
    cutoff = find(points > timebefore & points < (last-buffer));% remove = because ran into a problem when = timebefore
    point_cutoff = points(cutoff) + buffer; %throw away stimulus that is too short from the timebefore (added 20 to account for delay in the stimulus 20180125 figure out what this should be)
    ave = mean(stimulus(n,:));
    spiketriggered = zeros(1, timebefore + 1);
    if isempty(cutoff) == 0
        for cnt = 1:length(point_cutoff);
            spiketriggered = spiketriggered + ((stimulus(n, (point_cutoff(cnt)-timebefore):point_cutoff(cnt)))); % subtract off mean of stimulus to zero
        end
        STA(n, :) = (spiketriggered)./cnt;
    end
end

STA_mean = mean(STA);

figure
plot(STA_mean, 'k')

%% Flip STA to make the linear filter

clear Linear_filter

for cnt = 1:length(STA_mean) - 1
    Linear_filter(cnt) = STA_mean(timebefore + 2 - cnt);
end

figure
plot(Linear_filter)


%% To calculate generator signal. linear prediction of the response

clear generator
for n = 1:num
    generator(n, :) = conv(Linear_filter, stimulus(n,:));
end

%figure
%plot(generator(1, :))

%% To cut out  extraneous parts of the generator potential: beginning and end

prepoint = 5000;

Generator_potential = [];
Stimulus = [];
for n = 1:num
    Generator_potential = cat(2, Generator_potential, generator(n, timebefore:timebefore+length(stimulus)-1));
    Stimulus = cat(2, Stimulus, stimulus(n, prepoint+1:(length(stimulus) - prepoint)));
end

%% Normalize the linear filter (injection)

Generator_potential_SD = std(Generator_potential);
Stimulus_SD = std(Stimulus);
Stimulus_mean = mean(Stimulus);

Linear_filter_norm = (Linear_filter/Generator_potential_SD)*(Stimulus_SD);


figure
plot(Linear_filter_norm)

%% save and normalize

SampleInterval = 0.0001;

LN = Linear_filter_norm'./ (SampleInterval); % to normalize the filter by the sampling interval

figure
plot(LN)

%% Convolution again with the normalized linear filter

clear generator_norm
for n = 1:num
    generator_norm(n, :) = conv(Linear_filter_norm, stimulus(n,:));
end

Generator_potential_norm = [];

for n = 1:num
    Generator_potential_norm = cat(1, Generator_potential_norm, generator_norm(n, timebefore:timebefore+length(stimulus)-1));
end

%% Make a continous spike rate


SampleRate = 10000;

for n = 1:num
    M= [1:length(stimulus(n,:))] - 0.5;% try to make new bins for the histogram
    N(n, :) = hist(spikes.sp{n}+20, M); % 
    Spike_hist(n, :) = N(n,:)*SampleRate; % to turn the spike rate into time (seconds)
end

figure
plot(Spike_hist(n, :))


%% Cut off time before

output = [];
input = [];

for n = 1:num
    clear Spike_Hz
    Spike_Hz = Spike_hist(n, timebefore:(end-1));
    output = cat(2, Spike_hist(n, timebefore:(end-1)), output);
    input = cat(2, Generator_potential_norm(n, 1:length(Spike_Hz)), input);
end

%GP_norm = Generator_potential_norm(1:length(Spike_Hz));

figure
plot(input, output, 'o')



%% sort and bin

[G, index] = sort(input);
S = output(index);

bin = 10000;% added a 0 for current injection

clear G_ave S_ave

for cnt = 1:(length(output)./bin)
    G_ave(cnt) = mean(G((cnt-1)*bin+1:cnt*bin));
    S_ave(cnt) = mean(S((cnt-1)*bin+1:cnt*bin));
end

figure
plot(G_ave, S_ave, 'ko')


%% if outward peak

[peak, I] = max(LN);
tpeak = I*SampleInterval;% decimated by 10. Results in msec

%% Fit nonlinearity with line and take parameters

%Error bar calculation
% Fitting
I = [min(G_ave):abs((min(G_ave)/10)):max(G_ave)];

% Line fit
clear fit fitcoef
coef = [1 0];
[linefitcoef, resid, jacobian] = nlinfit(G_ave', S_ave', 'linefit', coef);
ConfidenceInterval = nlparci(linefitcoef, resid, jacobian)
Linefit = linefit(linefitcoef, I');

figure
clf
plot(G_ave, S_ave, 'ko', I, Linefit, 'r')
linefitcoef

%% Fit nonlinearity with line and take parameters only after a threshold

%Error bar calculation
% Fitting
I = [min(G_ave):abs((min(G_ave)/10)):max(G_ave)];

index = find(G_ave>2.8e-10 & G_ave <5e-10);
G_short = G_ave(index);
S_short = S_ave(index);

%% Line fit
clear fit fitcoef
coef = [1 0];
[linefitcoef, resid, jacobian] = nlinfit(G_short', S_short', 'linefit', coef);
ConfidenceInterval = nlparci(linefitcoef, resid, jacobian)
Linefit = linefit(linefitcoef, I');

figure
clf
plot(G_ave, S_ave, 'ko', I, Linefit, 'r')
linefitcoef
%%

setting = 'CurrenttoSpikesCone'


 %%
[power_x, power_y] = PowerSpectrumFinder(LN', SampleInterval);

figure
loglog(power_x, power_y)

%% return to this later
   psdest = psd(spectrum.periodogram,power_y,'Fs',SampleInterval,'NFFT',length(power_y)); % periodogram; try with default hanning
   normcumsumpsd = cumsum(psdest.Data)./sum(psdest.Data);
   Ind = find(normcumsumpsd <=0.5,1,'last');
   fprintf('Median frequency is %2.3f Hz\n',psdest.Frequencies(Ind));
%%
powerdensity = cumsum(power_y);
figure
loglog(power_x, powerdensity)

% power_median = psdest.Frequencies(Ind) % return to later

%%
linear.(setting) = LN;
timetopeak.(setting) = tpeak;
zeropower.(setting) = (power_y(1));
[peak, I] = max(power_y); 
maxpower.(setting)= power_x(I);
%medianfrequency.(setting) = power_median; % return to this later


generator_input.(setting) = G_ave';
current.(setting) = S_ave';
linefitslope.(setting) = linefitcoef(1);
linefitintercept.(setting) = linefitcoef(2);
offset.(setting) = mean(G_ave);

%%

I = (min(CellParameters.NonlinearInput.(setting))):1e-11:(max(CellParameters.NonlinearInput.(setting)));

% 20220602 Error: 
% Error using ==> nlinfit at 117
% The model function 'sigmoidbase' was not found.

% Sigmoid fit
clear fit sigfitcoef
coef = [0.2 200 5e-10 9e-11];
%coef = [0.2 250 5e-10 9e-11];
[sigfitcoef, resid, jacobian] = nlinfit((CellParameters.NonlinearInput.(setting)), (CellParameters.NonlinearOutput.(setting)), 'sigmoidbase', coef);
ConfidenceInterval = nlparci(real(sigfitcoef), real(resid), real(jacobian))
Sigfit = sigmoidbase(real(sigfitcoef), I');
 
hold on
plot((I), Sigfit,'r')

CellParameters.nonlinear_sigmoid.(setting) = real(sigfitcoef); % base, max, I1/2,  exponent


%% can export nonlinearity to Igor if you want
data = cat(2, CellParameters.NonlinearInput.(setting), CellParameters.NonlinearOutput.(setting));
