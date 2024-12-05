% Written by Scott Harris

%% Define the Direction Between Control and DTR

% These data are fit to differences of two gaussians over the same domain. 
% Then we convert these to these equivalent n-dimensional vectors. 
% The mean vector for both control and DTR is computed, and the difference between these mean vectors (which is itself n-dimensional) defines the direction that one must travel to get from a control filter to a DTR filter. 
% By projecting individual spatial filters onto this "difference vector," we can measure how similar they are to control versus DTR filters.

% The gaussian function is described by
% g(x) = ae^(-(x-2b)^2)/2c^2)

% Where a is the maximum value, b is the center location of the peak, and c is the standard deviation.
% Spatial filters are defined by a difference of two gaussian functions, one for the center and one for the surround:
% f(x) = gcenter(x) - gsurround(x)
% For practical purposes, we'll calculate the filter value discretely across a definite interval of [-1.52, 1.52] with 1000 observations. 
% Thus, in vector notion, a spatial filter is described by, where the subscript i denotes the cell.

%SET THE FILE NAMES THAT YOU WANT TO USE FOR EACH GROUP HERE
fitParamFile_cntrl = 'filename'; 
fitParamFile_exp = 'filename';

cntrlParams = readmatrix(fullfile(filePath, fitParamFile_cntrl));
expParams = readmatrix(fullfile(filePath, fitParamFile_exp));


%define the domain that you'd like to operate over
domain = linspace(-1.52, 1.52, 1000);

%define the gaussian function
gaussian = @(x, a, b, c) a*exp(-((x-b).^2)./(2*c^2));

% extract the mean parameters for each group. We'll use C for data for the =
% center gaussian and S for data for the surround gaussian

%compute the difference of gaussians for control and exp
controlFilter = gaussian(domain, mean(cntrlParams(:, 3)), 0, mean(cntrlParams(:, 1))) - gaussian(domain, mean(cntrlParams(:, 4)), 0, mean(cntrlParams(:, 2)));
experimentalFilter = gaussian(domain, mean(expParams(:, 3)), 0, mean(expParams(:, 1))) - gaussian(domain, mean(expParams(:, 4)), 0, mean(expParams(:, 2)));

differenceVector = experimentalFilter - controlFilter;

%first for control cells
cntrlDotProducts = zeros(1, size(cntrlParams, 1));
for i = 1:size(cntrlParams, 1)
    filter_i = gaussian(domain, cntrlParams(i, 3), 0, cntrlParams(i, 1)) - gaussian(domain, cntrlParams(i, 4), 0, cntrlParams(i, 2));
    differenceVector_i = filter_i - controlFilter;
    cntrlDotProducts(i) = dot(differenceVector_i, differenceVector);
end

%then repeat for experimental cells
expDotProducts = zeros(1, size(expParams, 1));
for i = 1:size(expParams, 1)
    filter_i = gaussian(domain, expParams(i, 3), 0, expParams(i, 1)) - gaussian(domain, expParams(i, 4), 0, expParams(i, 2));
    differenceVector_i = filter_i - controlFilter;
    expDotProducts(i) = dot(differenceVector_i, differenceVector);
end

%divide by norm squared, which gives you the projection normalized to the length of the difference vector, so that the difference vector itself should have a value of 1
normedProjection = @(x) x/norm(differenceVector)^2;


allVals = normedProjection([cntrlDotProducts, expDotProducts]); 
boxPlotGroups = zeros(1, numel(allVals));
boxPlotGroups(numel(cntrlDotProducts) + 1: end) = 1;

%% Shift the Control Filters To Look Like DTR and Recompute the Projection

% Now that we have a metric for determining how similar any given filter is to control and DTR, we will ask what the most important parameters are that distinguish the two conditions.
% We have chosen to examine the following parameters:
    % CA - center amplitude
    % SA - surround amplitude
    % CW - center width (standard deviation)
    % SW - surround width (standard deviation)
% Here, we take the control filters and set one parameter at a time to be equal to the mean value from DTR filters. We then ask how far along the difference vector this new "shifted filter" falls. In other words, how much does changing the parameter make the control filter look like a DTR filter?

maxC_shiftedSFs = zeros(size(cntrlParams, 1), size(controlFilter, 2));
maxS_shiftedSFs = zeros(size(cntrlParams, 1), size(controlFilter, 2));
wC_shiftedSFs = zeros(size(cntrlParams, 1), size(controlFilter, 2));
wS_shiftedSFs = zeros(size(cntrlParams, 1), size(controlFilter, 2));
allParams_shiftedSF = zeros(size(cntrlParams, 1), size(controlFilter, 2));

%The order in the csv is stdC, stdS, maxC, maxS
maxC_exp = mean(expParams(:, 3));
maxS_exp = mean(expParams(:, 4));
wC_exp = mean(expParams(:, 1));
wS_exp = mean(expParams(:, 2));

for i = 1:size(cntrlParams, 1)  % changing one exp parameter each for control SF 
    maxC_shiftedSFs(i, :) = gaussian(domain, maxC_exp, 0, cntrlParams(i, 1)) - gaussian(domain, cntrlParams(i, 4), 0, cntrlParams(i, 2));
    maxS_shiftedSFs(i, :) = gaussian(domain, cntrlParams(i, 3), 0, cntrlParams(i, 1)) - gaussian(domain, maxS_exp, 0, cntrlParams(i, 2));
    wC_shiftedSFs(i, :) = gaussian(domain, cntrlParams(i, 3), 0, wC_exp) - gaussian(domain, cntrlParams(i, 4), 0, cntrlParams(i, 2));
    wS_shiftedSFs(i, :) = gaussian(domain, cntrlParams(i, 3), 0, cntrlParams(i, 1)) - gaussian(domain, cntrlParams(i, 4), 0, wS_exp);
    allParams_shiftedSF(i, :) = gaussian(domain, maxC_exp, 0, wC_exp) - gaussian(domain, maxS_exp, 0, wS_exp);
end

%calculate difference vectors
maxC_differences = maxC_shiftedSFs - controlFilter;
maxS_differences = maxS_shiftedSFs - controlFilter;
wC_differences = wC_shiftedSFs - controlFilter;
wS_differences = wS_shiftedSFs - controlFilter;
allParams_differences = allParams_shiftedSF - controlFilter;

maxC_dotProducts = zeros(1, size(maxC_differences, 1));
maxS_dotProducts = zeros(size(maxC_dotProducts));
wC_dotProducts = zeros(size(maxC_dotProducts));
wS_dotProducts = zeros(size(maxC_dotProducts));
allParams_dotProducts = zeros(size(maxC_dotProducts));

for i = 1:numel(maxC_dotProducts)
    maxC_dotProducts(i) = dot(maxC_differences(i, :), differenceVector);
    maxS_dotProducts(i) = dot(maxS_differences(i, :), differenceVector);
    wC_dotProducts(i) = dot(wC_differences(i, :), differenceVector);
    wS_dotProducts(i) = dot(wS_differences(i, :), differenceVector);
    allParams_dotProducts(i) = dot(allParams_differences(i, :), differenceVector);
end

numControls = numel(cntrlDotProducts);
allVals = normedProjection([cntrlDotProducts, maxC_dotProducts, maxS_dotProducts, wC_dotProducts, wS_dotProducts, allParams_dotProducts, expDotProducts]);
boxPlotGroups = zeros(1, numel(allVals));
boxPlotGroups(numControls + 1: numControls*2) = 1;
boxPlotGroups(numControls*2 + 1: numControls*3) = 2;
boxPlotGroups(numControls*3 + 1: numControls*4) = 3;
boxPlotGroups(numControls*4 + 1: numControls*5) = 4;
boxPlotGroups(numControls*5 + 1: numControls*6) = 5;
boxPlotGroups(numControls*6 + 1:end) = 6;

%% Test All Combinations of Parameter Changes

% Next, we will test every combination of parameters, instead of just changing one at a time.
% This plot shows the amount that each combination of parameters moved control toward DTR. 
% Each point shows the average value across all control cells (i.e., just showing the mean here, not the full distribution).

aramNames = {'CA', 'SA', 'CW', 'SW'};
numParams = numel(paramNames);
totalCombos = 2^numel(paramNames);

expParams = [maxC_exp, maxS_exp, wC_exp, wS_exp];

% Preallocate matrix for combinations
allCombos = zeros(totalCombos, numParams);

ComboIDs = cell(totalCombos, 1);
meanProjectionPerCombo = zeros(totalCombos, 1);
% Loop through each combination
for c = 0:(totalCombos-1)
    % Convert index c to binary vector to determine which elements to replace
    binVec = fliplr(dec2bin(c, numParams) - '0'); % '0' and '1' -> 0 or 1
    % Start with defaultValues
    row = [0 0 0 0];
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
    projections_c = zeros(size(cntrlParams, 1), 1);
    for i = 1:size(cntrlParams, 1)

        params_i = cntrlParams(i, [3 4 1 2]);
        for j = 1:numel(row)
            if row(j)
                params_i(j) = expParams(j);
            end
        end

        shiftedFilter = gaussian(domain, params_i(1), 0, params_i(3)) - gaussian(domain, params_i(2), 0, params_i(4));
        projections_c(i) = normedProjection(dot(shiftedFilter - controlFilter, differenceVector));
    end
    meanProjectionPerCombo(c+1) = mean(projections_c, 'omitnan');
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
