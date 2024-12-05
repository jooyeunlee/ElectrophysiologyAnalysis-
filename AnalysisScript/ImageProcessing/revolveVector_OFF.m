% Written by Luca Della Santina

function [RFout, RFoutSigned] = revolveVector(RFin)
% Move the peak of the RF in the middle
% [~,iMax] = max(RFin);

% 20230721 changed max to abs(min) because of OFF RF - to normalize to the
% negative peak
[~,iMin] = (min(RFin));
RFinShifted = zeros(numel(RFin),1);
if iMin > numel(RFin)/2
    % Signal neds to be shifted leftwards
    RFinShifted(1:ceil(numel(RFin)/2)) = RFin(iMin-floor(numel(RFin)/2):iMin);
    newMid = iMin-(iMin-floor(numel(RFin)/2))+2;
    RFinShifted(newMid:newMid+numel(RFin)-iMin-1) = RFin(iMin+1:numel(RFin));
    RFinShifted(numel(RFin)+ceil(numel(RFin)/2)-iMin+1:end) = RFin(1:iMin- ceil(numel(RFin)/2));
    RFin = RFinShifted;
end

% Make the RF symmetric by averaging the left and right sides
RFinNorm = RFin - max(RFin);
RFinNorm = RFinNorm/abs(min(RFinNorm));
RFinNormFlipped = fliplr(RFinNorm')';
RFsmooth = zeros(numel(RFinNorm),1);
for i = 1:numel(RFinNorm)
        RFsmooth(i) = (RFinNorm(i) + RFinNormFlipped(i)) /2;
end
RFsmoothNorm = RFsmooth - max(RFsmooth);
RFsmoothNorm = RFsmoothNorm/abs(min(RFsmoothNorm));

% plot(RFinNorm, '-k');
% hold on
% plot(RFinNormFlipped, '-b');
% plot(RFsmoothNorm, '-r');
% legend({'RFin'; 'RFin flipped'; 'RFin smoothed'});

%% Revolve the 1D vector to have a 2D image
v1D = RFsmoothNorm;
% create a meshgrid of n x n given a 1D vector of length n  
[rx,ry] = meshgrid( -floor(numel(v1D)/2):floor(numel(v1D)/2) , -floor(numel(v1D)/2):floor(numel(v1D)/2) );
% calculate the radius for each point of the meshgrid
r = hypot( rx , ry );
% extract half of the RFsmooth signal
halfRFsmooth = v1D((length(v1D)-1)/2 + 1 : end );
% interpolate the missing values for each radius distance
vq = interp1( 0:(length(halfRFsmooth)-1) , halfRFsmooth , r(:) , 'linear' , 0 );
% write results in RFout
RFout = r; 
RFout(:) = vq;
RFout = imresize(RFout,[67,67]);

% Zero the far periphery of the RF (initial values of the vector)
RFoutSigned = RFout - mean(RFsmoothNorm(1:5));
% Re-normalize between baseline and max
RFoutSigned = RFoutSigned./abs(min(RFoutSigned(:)));

% imshow(uint8(222*RFout/max(RFout(:))))
end