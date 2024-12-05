function [RFSTout] = revolveVector3(RFin, sz, ONcell)
%% Revolves a 1D spatial RF into a 2D spatial RF for each time of the spatio-temporal input RF
%
% Input:
% RFin = spatio-temporal nxm RF in which:
%        rows n are 1D spatial receptive fields 
%        columns m are each time point of the temporal RF
%
% sz [n,m] = optional size of the output spatial RF resized image,
%              if not specified, output size will be same as input
% Output:
% RFSTout = 3D n x m x t spatio-temporal RF with:
%           n x m = Spatial RF converted to 2D images by revolution
%           t = times of the temporal RF
%
% Spatial RF values are normalized to max of original spatial RF

% % Normalize the RF to peak value 
% if OFFcell
%     % OFF cells, Normalize the RF from baseline=0 to negative peak=-1
%     RFin(:) = RFin(:)/abs(min(RFin(:))); % Original OFF RF looks like ON cell and this is to flip polarity of OFF RF to move the pick below 0
%    
% else
%     % ON cells, Normalize the RF from baseline=0 to positive peak=1
%     RFin(:) = RFin(:)/max(RFin(:));
% end


% For analyzing ON cell 
% 20240515: The polarity of RF has been flipped as OFF cells in Igor, Data were saved as "ON_RF_Flipped_data.mat", Treat this the same as OFF RF

if ONcell
%    % ON cells, Normalize the RF from baseline=0 to positive peak=1
%    RFin(:) = RFin(:)/max(RFin(:));
    RFin(:) = RFin(:)/abs(min(RFin(:)));
else
%    % OFF cells
%    RFin(:) = RFin(:)/abs(min(RFin(:)));
    RFin(:) = RFin(:)/max(RFin(:));
end


for t = 1:size(RFin,2)

%     % Normalize the RF to peak value
%     if OFFcell
%         % OFF cells, Normalize the RF from baseline=0 to negative peak=-1
%         RFin(:,t) = RFin(:,t)/abs(min(RFin(:,t)));
%     else
%         % ON cells, Normalize the RF from baseline=0 to positive peak=1
%         RFin(:,t) = RFin(:,t)/max(RFin(:,t));        
%     end

    % Revolve the 1D vector to have a 2D image
    v1D = RFin(:,t);

    % create a meshgrid of n x n given a 1D vector of length n
    [rx,ry] = meshgrid( -floor(numel(v1D)/2):floor(numel(v1D)/2) , -floor(numel(v1D)/2):floor(numel(v1D)/2) );

    % calculate the radius for each point of the meshgrid
    r = hypot( rx , ry );

    % extract half of the RFsmooth signal
    % halfRFsmooth = v1D((length(v1D)-1)/2 + 1 : end );
    % interpolate the missing values for each radius distance
    vq = interp1( 0:(length(v1D)-1) , v1D , r(:) , 'linear' , 0 );
    % write results in RFout
    RFout = r;
    RFout(:) = vq;
    
    % Resize the receptive field instead of decimating it
    % Use this if you want the output RF be a different size than input
    % RFout = imresize(RFout,[67,67]);
    if ~isempty(sz)
        RFout = imresize(RFout,sz);
    end

    RFSTout(:,:,t) = RFout;

%     % Show individual 2D spatial RF image at time t
%     imshow(uint8(255*RFout/max(RFout(:))))
end

% % Show final RFST movie
% for i=1:size(RFSTout,3)
%     imshow(RFSTout(:,:,i));
% end
end