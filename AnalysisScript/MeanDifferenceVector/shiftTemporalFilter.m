% Written by Scott Harris

function [shiftedFilter] = shiftTemporalFilter(temporalFilter, timeToPeakShift_ms, timeToTroughShift_ms, peakAmplitudeScale, troughAmplitudeScale, zeroCrossingShift_ms, varargin)
%shiftTemporalFilter Shifts the temporal filter of a cell according to the
%input paramets
%   INPUTS:
% - filter: timeseries of the original temporal filter
% - timeToPeakShift_ms: number of milliseconds to shift the peak by
% - timeToTroughShift_ms: number of milliseconds to shift the trough by
% - peakAmplitudeChange: scale factor for the peak amplitude. E.g., 1 means keep the original value, 2 means double it, -1 means invert it
% - troughAmplitudeChange: scale factor for the trough amplitude. E.g., 1 means keep the originalvalue, 2 means double it, -1 means invert it
% - zeroCrossingShift_ms: number of milliseconds to shift the zero crossing time by
% - varargin:
%     - "sampleRate": in hertz, the sample rate of the temporal filter. Default is 30
%     - "useAbsoluteTimes": boolean 1 or 0. If 1, then
% timeToPeakShift_ms and timeToTroughShift_ms are given as the
% actual time points of the new peak/trough, not as deltas. Units
% are assumed to be in seconds in this case
%   OUTPUTS:
% - shiftedFilter: the temporal filter shifted in time and size based on the inputs

p = inputParser;
defaultSampleRate = 100; % hz
addParameter(p,'sampleRate',defaultSampleRate,@isnumeric)
p.CaseSensitive = 0;
addParameter(p, 'useAbsoluteTimes', 0, @islogical)
addParameter(p, 'makeImage', 0, @islogical)
parse(p, varargin{:});

sampleRate = p.Results.sampleRate;
useAbsoluteTimes = p.Results.useAbsoluteTimes;
makeImage = p.Results.makeImage;

dwell = 1/sampleRate; %seconds per frame

timestamps = 0:dwell:dwell*(numel(temporalFilter)-1); %get the timestamp for each value in the filter

[originalPeak, originalPeakIndex] = max(temporalFilter);
[originalTrough, originalTroughIndex] = min(temporalFilter);
originalTTP = timestamps(originalPeakIndex);
originalTTT = timestamps(originalTroughIndex);

%convert to seconds
TTPShift_seconds = timeToPeakShift_ms/1000;
TTTShift_seconds = timeToTroughShift_ms/1000;
TTZShift_seconds = zeroCrossingShift_ms/1000;

%calculate existing zero crossing time (must occur after the peak and
%before the trough)
btwnTroughs = temporalFilter(originalPeakIndex:originalTroughIndex);
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


if useAbsoluteTimes
    TTPShift_seconds = timeToPeakShift_ms - originalTTP; %in this case, timeToPeakShift_ms is assumed to actually be in seconds, not milliseconds
    TTTShift_seconds = timeToTroughShift_ms - originalTTT;
    TTZShift_seconds = zeroCrossingShift_ms - originalZeroCrossing;
    disp('absolute values')
end


%first handle the amplitude changes by scaling each value in the filter
%based on a linear combination of how close it is to the peak or trough
distanceToPeak = (originalPeak - temporalFilter)./(originalPeak - originalTrough); % the smaller this number is to 0 the closer the value is to the peak. The closer it is to 1, the closer it is to the trough
rescaledFilter = temporalFilter.*(((1 - distanceToPeak) * peakAmplitudeScale) + distanceToPeak * troughAmplitudeScale);
%for the tail bit, scale based on a linear combination of how close each point is to
%the trough vs ending point timestamp
rescaledEndSegment = rescaledFilter(originalTroughIndex+1:end);
originalEndSegment = temporalFilter(originalTroughIndex+1:end);
proximityToEnd = (timestamps(originalTroughIndex+1:end)-timestamps(originalTroughIndex))/(timestamps(end)-timestamps(originalTroughIndex));
mixedEndSegment = rescaledEndSegment'.*(1-proximityToEnd) + originalEndSegment' .* proximityToEnd;
rescaledFilter(originalTroughIndex+1:end) = mixedEndSegment;

%next shift the peak, zero crossing, and trough time:

%break the trace up into 4 segments:
xi1 = timestamps(1:originalPeakIndex);
xi2 = timestamps(originalPeakIndex+1:originalZeroCrossingIndex);
xi3 = timestamps(originalZeroCrossingIndex+1:originalTroughIndex);
xi4 = timestamps(originalTroughIndex+1:end);


w1 = 1-(originalTTP - xi1)/originalTTP;
xf1 = xi1 + w1 * TTPShift_seconds;

w2 = 1-(originalTTZ-xi2)/(originalTTZ - originalTTP);
xf2 = xi2 + (1-w2)*TTPShift_seconds + w2 * TTZShift_seconds;

w3 = 1-(originalTTT - xi3)/(originalTTT - originalTTZ);
xf3 = xi3 + (1-w3)*TTZShift_seconds + w3*TTTShift_seconds;

w4 = 1 - (timestamps(end) - xi4)/(timestamps(end) - originalTTT);
xf4 = xi4 + (1-w4)*TTTShift_seconds;

newTimestamps = [xf1, xf2, xf3, xf4];

if ~sum(sort(newTimestamps) == newTimestamps) == numel(newTimestamps)
    disp('Could not shift the filter due to incompatible requirements')
    shiftedFilter = NaN;
    return
end

shiftedFilter = interp1(newTimestamps, rescaledFilter, timestamps); %recreate a temporal filter that can be used with the original timestamps

if makeImage
    figure
    title('Shifted Temporal Filter')
    hold on
    plot(timestamps, temporalFilter, '-k', 'LineWidth', 1.5)
    plot(timestamps, shiftedFilter, '-r', 'LineWidth', 1.5)
    xlabel('Time')
    ylabel('Amplitude')
    l = legend({'Original Filter', 'New Filter'});
    l.Location = 'best';
end

end