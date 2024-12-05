% Written by Scott Harris

%%  Define the Direction Between Control and DTR

% Each temporal filter is a time series of 50 points. We turn these time series into equivalent 50-dimensional vectors (or points in 50-dimensional space). 
% The mean vector for both control and DTR is then computed. The difference between these mean vectors (which is itself a 50-dimensional vector) defines the direction that one must travel to get from a control filter to a DTR filter. 
% By projecting individual temporal filters onto this "difference vector," we can measure how similar they are to control versus DTR filters.

% Set the paramters 
% TTP_=;
% TTT_=;
% PeakAmp_=;
% TroughAmp_=;
% TTZ=;

sampleRate = 100; %sample rate of temporal filters

% reorganize the data
numCells1 = size(data1, 1)/50;
controlTFs = zeros(50, numCells1);
for i = 1:numCells1
    TF_i = data1(i:numCells1:size(data1, 1));
    controlTFs(:, i) = TF_i;
end

numCells2 = size(data2, 1)/50;
expTFs = zeros(50, numCells2);
for i = 1:numCells2
    TF_i = data2(i:numCells2:size(data2, 1));
    expTFs(:, i) = TF_i;
end

meanExpFilter = mean(expTFs, 2);
meanControlFilter = mean(controlTFs, 2);
meanDifferenceVector = meanExpFilter - meanControlFilter;

%calculate the dot product between each filter and the difference filter
controlDifferenceVectors = controlTFs - meanControlFilter;
expDifferenceVectors = expTFs - meanControlFilter;

controlDotProducts = zeros(1, size(controlDifferenceVectors, 2));
expDotProducts = zeros(1, size(expDifferenceVectors, 2));

for i = 1:numel(controlDotProducts)
    controlDotProducts(i) = dot(controlDifferenceVectors(:, i), meanDifferenceVector);
end
for i = 1:numel(expDotProducts)
    expDotProducts(i) = dot(expDifferenceVectors(:, i), meanDifferenceVector);
end

normedProjection = @(x) x./norm(meanDifferenceVector)^2;

allVals = normedProjection([controlDotProducts, expDotProducts]);
boxPlotGroups = zeros(1, numel(allVals));
boxPlotGroups(numel(controlDotProducts) + 1: end) = 1;

p = ranksum(controlDotProducts, expDotProducts)

%% Shift the Control Filters To Look Like DTR and Recompute the Projection

% Now that we have a metric for determining how similar any given filter is to control and DTR, we will ask what the most important parameters are that distinguish the two conditions.
% We have chosen to examine the following parameters:
    % TTP - time to peak
    % TTT - time to trough
    % PA - peak amplitude
    % TA - trough amplitude
    % TTZ - time to zero crossing
% Here, we take the control filters and set one parameter at a time to be equal to the mean value from DTR filters. We then ask how far along the difference vector this new "shifted filter" falls. 
% In other words, how much does changing the parameter make the control filter look like a DTR filter?

ttp_shiftedTFs = zeros(size(controlTFs));
ttt_shiftedTFs = zeros(size(controlTFs));
peakAmp_shiftedTFs = zeros(size(controlTFs));
troughAmp_shiftedTFs = zeros(size(controlTFs));
ttz_shiftedTFs = zeros(size(controlTFs));

allParams_shiftedTFs = zeros(size(controlTFs));

dwell = 1/sampleRate; %seconds per frame
timestamps = 0:dwell:dwell*(size(controlTFs,1)-1); %get the timestamp for each value in the filter
for i = 1:size(controlTFs, 2)
    tf_i = controlTFs(:, i);
    
    [originalPeak, originalPeakIndex] = max(tf_i);
    [originalTrough, originalTroughIndex] = min(tf_i);
    originalTTP = timestamps(originalPeakIndex);
    originalTTT = timestamps(originalTroughIndex);
    
    btwnTroughs = tf_i(originalPeakIndex:originalTroughIndex);
    timestampsBtwnTroughs = timestamps(originalPeakIndex:originalTroughIndex);
    [~, nearZeroIndex] = min(abs(btwnTroughs));
    originalZeroCrossingIndex = nearZeroIndex + originalPeakIndex - 1;

    v1 = btwnTroughs(nearZeroIndex);
    t1 = timestampsBtwnTroughs(nearZeroIndex);
    if v1 > 0
        v2 = btwnTroughs(nearZeroIndex + 1);
        t2 = timestampsBtwnTroughs(nearZeroIndex + 1);
    else
        v2 = btwnTroughs(nearZeroIndex - 1);
        t2 = timestampsBtwnTroughs(nearZeroIndex - 1);
    end
    originalTTZ = (1-abs(v1/(v2-v1)))*t1 + abs(v1/(v2-v1))*t2;

    TTPShift = (TTP_exp - originalTTP)*1000;
    TTTShift = (TTT_exp - originalTTT)*1000;
    TTZShift = (TTZ_exp - originalTTZ)*1000;
    peakAmpShift = PeakAmp_exp/originalPeak;
    troughAmpShift = TroughAmp_exp/originalTrough;

    ttp_shiftedTFs(:, i) = shiftTemporalFilter(tf_i, TTPShift, 0, 1, 1, 0, 'sampleRate', sampleRate, 'makeImage', false);
    ttt_shiftedTFs(:, i) = shiftTemporalFilter(tf_i, 0, TTTShift, 1, 1, 0, 'sampleRate', sampleRate, 'makeImage', false);
    peakAmp_shiftedTFs(:, i) = shiftTemporalFilter(tf_i, 0, 0, peakAmpShift, 1, 0, 'sampleRate', sampleRate, 'makeImage', false);
    troughAmp_shiftedTFs(:, i) = shiftTemporalFilter(tf_i, 0, 0, 1, troughAmpShift, 0, 'sampleRate', sampleRate, 'makeImage', false);
    ttz_shiftedTFs(:, i) = shiftTemporalFilter(tf_i, 0, 0, 1, 1, TTZShift, 'sampleRate', sampleRate', 'makeImage', false);
    allParams_shiftedTFs(:, i) = shiftTemporalFilter(tf_i, TTPShift, TTTShift, peakAmpShift, troughAmpShift, TTZShift, 'sampleRate', sampleRate, 'makeImage', false);
end

%calculate the difference vectors
ttpDifferenceVectors = ttp_shiftedTFs - meanControlFilter;
tttDifferenceVectors = ttt_shiftedTFs - meanControlFilter;
peakAmpDifferenceVectors = peakAmp_shiftedTFs - meanControlFilter;
troughAmpDifferenceVectors = troughAmp_shiftedTFs - meanControlFilter;
ttzDifferenceVectors = ttz_shiftedTFs - meanControlFilter;
allParamsDifferenceVectors = allParams_shiftedTFs - meanControlFilter;

ttpDotProducts = zeros(1, size(controlDifferenceVectors, 2));
tttDotProducts = zeros(1, size(tttDifferenceVectors, 2));
peakAmpDotProducts = zeros(1, size(peakAmpDifferenceVectors, 2));
troughAmpDotProducts = zeros(1, size(troughAmpDifferenceVectors, 2));
ttzDotProducts = zeros(1, size(ttzDifferenceVectors, 2));
allParamsDotProducts = zeros(1, size(allParamsDifferenceVectors, 2));

for i = 1:numel(ttpDotProducts)
    ttpDotProducts(i) = dot(ttpDifferenceVectors(:, i), meanDifferenceVector);
end
for i = 1:numel(tttDotProducts)
    tttDotProducts(i) = dot(tttDifferenceVectors(:, i), meanDifferenceVector);
end
for i = 1:numel(peakAmpDotProducts)
    peakAmpDotProducts(i) = dot(peakAmpDifferenceVectors(:, i), meanDifferenceVector);
end
for i = 1:numel(troughAmpDotProducts)
    troughAmpDotProducts(i) = dot(troughAmpDifferenceVectors(:, i), meanDifferenceVector);
end
for i = 1:numel(ttzDotProducts)
    ttzDotProducts(i) = dot(ttzDifferenceVectors(:, i), meanDifferenceVector);
end
for i = 1:numel(allParamsDotProducts)
    allParamsDotProducts(i) = dot(allParamsDifferenceVectors(:, i), meanDifferenceVector);
end

numControls = numel(controlDotProducts);
allVals = normedProjection([controlDotProducts, ttpDotProducts, tttDotProducts, peakAmpDotProducts, troughAmpDotProducts, ttzDotProducts, allParamsDotProducts, expDotProducts]);
boxPlotGroups = zeros(1, numel(allVals));
boxPlotGroups(numControls + 1: numControls*2) = 1;
boxPlotGroups(numControls*2 + 1: numControls*3) = 2;
boxPlotGroups(numControls*3 + 1: numControls*4) = 3;
boxPlotGroups(numControls*4 + 1: numControls*5) = 4;
boxPlotGroups(numControls*5 + 1: numControls*6) = 5;
boxPlotGroups(numControls*6 + 1: numControls*7) = 6;
boxPlotGroups(numControls*7 + 1:end) = 7;

p_ControlVsAllParams = ranksum(controlDotProducts, allParamsDotProducts)

%% Test All Combinations of Parameter Changes

% Next, we will test every combination of parameters, instead of just changing one at a time.
% This plot shows the amount that each combination of parameters moved control toward DTR. Each point shows the average value across all control cells (i.e., just showing the mean here, not the full distribution).

defaultValues = [0, 0, 1, 1, 0]; %default values to enter into the shift function
paramNames = {'TTP', 'TTT', 'PA', 'TA', 'TTZ'};
numParams = numel(defaultValues);
totalCombos = 2^numel(defaultValues);

% Preallocate matrix for combinations
allCombos = zeros(totalCombos, numParams);

ComboIDs = cell(totalCombos, 1);
meanProjectionPerCombo = zeros(totalCombos, 1);
% Loop through each combination
for c = 0:(totalCombos-1)
    % Convert index c to binary vector to determine which elements to replace
    binVec = fliplr(dec2bin(c, numParams) - '0'); % '0' and '1' -> 0 or 1
    % Start with defaultValues
    row = [0 0 0 0 0];
    id = '';
    
    % Replace elements in row based on binVec
    for j = 1:numParams
        if binVec(j) == 1
            row(j) = 1; 
            id = [id '-' paramNames{j}];
        end
    end
   
    allCombos(c+1, :) = row;
    ComboIDs{c+1} = id;

    %compute the shifted filters
    projection_c = zeros(size(controlTFs, 2), 1);
    for i = 1:size(controlTFs, 2)
        tf_i = controlTFs(:, i);

        [originalPeak, originalPeakIndex] = max(tf_i);
        [originalTrough, originalTroughIndex] = min(tf_i);
        originalTTP = timestamps(originalPeakIndex);
        originalTTT = timestamps(originalTroughIndex);

        btwnTroughs = tf_i(originalPeakIndex:originalTroughIndex);
        timestampsBtwnTroughs = timestamps(originalPeakIndex:originalTroughIndex);
        [~, nearZeroIndex] = min(abs(btwnTroughs));
        originalZeroCrossingIndex = nearZeroIndex + originalPeakIndex - 1;

        v1 = btwnTroughs(nearZeroIndex);
        t1 = timestampsBtwnTroughs(nearZeroIndex);
        if v1 > 0
            v2 = btwnTroughs(nearZeroIndex + 1);
            t2 = timestampsBtwnTroughs(nearZeroIndex + 1);
        else
            v2 = btwnTroughs(nearZeroIndex - 1);
            t2 = timestampsBtwnTroughs(nearZeroIndex - 1);
        end
        originalTTZ = (1-abs(v1/(v2-v1)))*t1 + abs(v1/(v2-v1))*t2;

        TTPShift = (TTP_exp - originalTTP)*1000;
        TTTShift = (TTT_exp - originalTTT)*1000;
        TTZShift = (TTZ_exp - originalTTZ)*1000;
        peakAmpShift = PeakAmp_exp/originalPeak;
        troughAmpShift = TroughAmp_exp/originalTrough;

        expParams = [TTPShift, TTTShift, peakAmpShift, troughAmpShift, TTZShift];
        params_i = [0 0 0 0 0];
        for j = 1:numel(row)
            expValues = [TTP_exp, TTT_exp, PeakAmp_exp, TroughAmp_exp, TTZ_exp];
            if row(j)
                params_i(j) = expParams(j);
            else
                params_i(j) = defaultValues(j);
            end
        end

        shiftedFilter = shiftTemporalFilter(tf_i, params_i(1), params_i(2), params_i(3), params_i(4), params_i(5), 'sampleRate', sampleRate, 'makeImage', false);
        projection_c(i) = normedProjection(dot(shiftedFilter' - meanControlFilter, meanDifferenceVector));
    end
    meanProjectionPerCombo(c+1) = mean(projection_c, 'omitnan');
end
ComboIDs{1} = 'None';


% calculate the marginal benefit of each parameter
meanMarginalAdditiveBenefits = zeros(1, numParams);
for i = 1:numParams
    additiveBenefits = [];
    for j = 1:totalCombos
        combo_j = allCombos(j, :);
        if ~combo_j(i)
            continue
        end
        dp_j = meanProjectionPerCombo(j);
        
        combo_noI = combo_j;
        combo_noI(i) = 0;
        row_noI = find(ismember(allCombos, combo_noI, 'rows'));
        
        dp_noI = meanProjectionPerCombo(row_noI);

        additiveBenefits(end+1) = dp_j - dp_noI;
    end
    meanMarginalAdditiveBenefits(i) = mean(additiveBenefits);
end

