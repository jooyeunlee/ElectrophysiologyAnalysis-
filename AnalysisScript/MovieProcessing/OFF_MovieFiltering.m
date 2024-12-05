%% Spatio-temporal Receptive Field filtering of an arbitrary B/W videoclip
%
% Version 1.0 2023-05-26 by Luca Della Santina
%
% + Load JL data OFF-S, OFF-T RGCs in Full, Partial, DTR conditions
% + Calculates the average stRF normalized by the negative peak response
% + Applies the stRF as kernel to convolve an arbirary movie
% + Calculates SSIM metrics for the convolved vs original videos
%
% Version 1.1 2023-05-28 by Luca Della Santina
% 
% + Added MSE calculation
% + Added PSNR calculation
% - Imposed conversion to 8bit of the original and processed video
% - Restricted SSIM calculation to ignore the last frames of the video
%     by [size of the temporal kernel] frames, to avoid polarity shift
%     due to last part of video filtered only with RF periphery

% Create lists for all the available data
lstFull_S = {'Full_S_RF0', 'Full_S_RF1', 'Full_S_RF2', 'Full_S_RF3', 'Full_S_RF4', 'Full_S_RF5', 'Full_S_RF6', 'Full_S_RF7', 'Full_S_RF8', 'Full_S_RF9', 'Full_S_RF10', 'Full_S_RF11', 'Full_S_RF12', 'Full_S_RF13', 'Full_S_RF14', 'Full_S_RF15', 'Full_S_RF16'};
lstFull_T = {'Full_T_RF0', 'Full_T_RF1', 'Full_T_RF2', 'Full_T_RF3', 'Full_T_RF4', 'Full_T_RF5', 'Full_T_RF6', 'Full_T_RF7', 'Full_T_RF8', 'Full_T_RF9', 'Full_T_RF10', 'Full_T_RF11', 'Full_T_RF12'};
lstDTR_S = {'DTR_S_RF0', 'DTR_S_RF1', 'DTR_S_RF2', 'DTR_S_RF3', 'DTR_S_RF4', 'DTR_S_RF5', 'DTR_S_RF6', 'DTR_S_RF7', 'DTR_S_RF8', 'DTR_S_RF9', 'DTR_S_RF10', 'DTR_S_RF11', 'DTR_S_RF12', 'DTR_S_RF13', 'DTR_S_RF14', 'DTR_S_RF15'};
lstDTR_T = {'DTR_T_RF0', 'DTR_T_RF1', 'DTR_T_RF2', 'DTR_T_RF3', 'DTR_T_RF4', 'DTR_T_RF5', 'DTR_T_RF6', 'DTR_T_RF7', 'DTR_T_RF8', 'DTR_T_RF9', 'DTR_T_RF10', 'DTR_T_RF11', 'DTR_T_RF12', 'DTR_T_RF13', 'DTR_T_RF14', 'DTR_T_RF15'};
lstPartial_S = {'Partial_S_RF0', 'Partial_S_RF1', 'Partial_S_RF2', 'Partial_S_RF3', 'Partial_S_RF4', 'Partial_S_RF5', 'Partial_S_RF6', 'Partial_S_RF7', 'Partial_S_RF8', 'Partial_S_RF9', 'Partial_S_RF10', 'Partial_S_RF11', 'Partial_S_RF12', 'Partial_S_RF13', 'Partial_S_RF14', 'Partial_S_RF15', 'Partial_S_RF16'};
lstPartial_T = {'Partial_T_RF0', 'Partial_T_RF1', 'Partial_T_RF2', 'Partial_T_RF3', 'Partial_T_RF4', 'Partial_T_RF5', 'Partial_T_RF6', 'Partial_T_RF7', 'Partial_T_RF8', 'Partial_T_RF9', 'Partial_T_RF10', 'Partial_T_RF11', 'Partial_T_RF12'};

% Average individual RFs per condition (note need to transpose vectors
% because RFs are [Time x Space] while revolveVector needs [Space X Time]

% sz = [9,9]; % Resize all spatial filters, leave empty [] for original size
sz = [71,71]; % to make RF size 1/18 of 1280x720 image size 
OFFrgc = true;

tmpRF = [];
for i = 1:numel(lstFull_S)
    tmpRF(:,:,:,i) = revolveVector3(eval(lstFull_S{i})', sz, OFFrgc);
end
% Average RFs and normalize to minimum = -1 for OFF RGCs
RF_Full_S_Avg = mean(tmpRF,4);
RF_Full_S_Avg(:) = RF_Full_S_Avg(:)/abs(min(RF_Full_S_Avg(:)));
%RF_Full_S_Avg = flip(RF_Full_S_Avg,3);

tmpRF = [];
for i = 1:numel(lstFull_T)
    tmpRF(:,:,:,i) = revolveVector3(eval(lstFull_T{i})', sz, OFFrgc);
end
RF_Full_T_Avg = mean(tmpRF,4);
RF_Full_T_Avg(:) = RF_Full_T_Avg(:)/abs(min(RF_Full_T_Avg(:)));


tmpRF = [];
for i = 1:numel(lstDTR_S)
    tmpRF(:,:,:,i) = revolveVector3(eval(lstDTR_S{i})', sz, OFFrgc);
end
RF_DTR_S_Avg = mean(tmpRF,4);
RF_DTR_S_Avg(:) = RF_DTR_S_Avg(:)/abs(min(RF_DTR_S_Avg(:)));


tmpRF = [];
for i = 1:numel(lstDTR_T)
    tmpRF(:,:,:,i) = revolveVector3(eval(lstDTR_T{i})', sz, OFFrgc);
end
RF_DTR_T_Avg = mean(tmpRF,4);
RF_DTR_T_Avg(:) = RF_DTR_T_Avg(:)/abs(min(RF_DTR_T_Avg(:)));


tmpRF = [];
for i = 1:numel(lstPartial_S)
    tmpRF(:,:,:,i) = revolveVector3(eval(lstPartial_S{i})', sz, OFFrgc);
end
RF_Partial_S_Avg = mean(tmpRF,4);
RF_Partial_S_Avg(:) = RF_Partial_S_Avg(:)/abs(min(RF_Partial_S_Avg(:)));


tmpRF = [];
for i = 1:numel(lstPartial_T)
    tmpRF(:,:,:,i) = revolveVector3(eval(lstPartial_T{i})', sz, OFFrgc);
end
RF_Partial_T_Avg = mean(tmpRF,4);
RF_Partial_T_Avg(:) = RF_Partial_T_Avg(:)/abs(min(RF_Partial_T_Avg(:)));


%% Show final sptatio-temporal RF as a movie
tmpRF = RF_Full_T_Avg;
disp(['Max value = ' num2str(max(tmpRF(:)))]);
disp(['Min value = ' num2str(min(tmpRF(:)))]);

for i=1:size(tmpRF,3)
    imshow(255*tmpRF(:,:,i),jet);
end

clear lst* tmp* i sz OFFrgc

%% Filter video file with average spatio-temporal RF for each condition

% Load the video file (requires Windows, Linux OS gives codecs error)
videofilename = ['1_720_1x.mp4'];

V = loadVideoFile(videofilename);

%% Batch analysis (Need analyzeBatchGroup.m and analyzeGroup.m)
Kernels = cat(4, RF_Full_S_Avg, RF_Partial_S_Avg, RF_DTR_S_Avg, RF_Full_T_Avg, RF_Partial_T_Avg, RF_DTR_T_Avg);
groupNames = {'Full_OFF_S', 'Partial_OFF_S', 'DTR_OFF_S', 'Full_OFF_T', 'Partial_OFF_T', 'DTR_OFF_T'};
[MSEBatch, PSNRBatch, SSIMBatch, diffErr] = analyzeBatchGroups(V, Kernels, videofilename, groupNames);

SSIMbatch_avg = mean(SSIMBatch,2)';
MSEbatch_avg  = mean(MSEBatch,2)';
PSNRbatch_avg = mean(PSNRBatch,2)';
diffErr_avg = mean(diffErr,2)'; 

