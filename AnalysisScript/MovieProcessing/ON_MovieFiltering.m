% Written by Luca Della Santina
% Edited by Jeremiah V John

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


% For ON cell
lstFull_S = {'Full_RF0', 'Full_RF1', 'Full_RF2', 'Full_RF3', 'Full_RF4', 'Full_RF5', 'Full_RF6', 'Full_RF7', 'Full_RF8', 'Full_RF9', 'Full_RF10', 'Full_RF11'};
lstDTR_S = {'DTR_RF0', 'DTR_RF1', 'DTR_RF2', 'DTR_RF3', 'DTR_RF4', 'DTR_RF5', 'DTR_RF6', 'DTR_RF7', 'DTR_RF8', 'DTR_RF9', 'DTR_RF10'};
lstPartial_S = {'Partial_RF0', 'Partial_RF1', 'Partial_RF2', 'Partial_RF3', 'Partial_RF4', 'Partial_RF5', 'Partial_RF6', 'Partial_RF7', 'Partial_RF8', 'Partial_RF9', 'Partial_RF10', 'Partial_RF11'};

% Average individual RFs per condition (note need to transpose vectors
% because RFs are [Time x Space] while revolveVector needs [Space X Time]
sz = [53,53]; % to make RF size 1/18 of 960x540 image size = 1/18.46
ONrgc = true;

tmpRF = [];
for i = 1:numel(lstFull_S)
    tmpRF(:,:,:,i) = revolveVector3(eval(lstFull_S{i})', sz, ONrgc);
end
% Average RFs and normalize to minimum = -1 for OFF RGCs (20240515: Treating ON cell
% the same as OFF cell as the polarity of ON RF has been flipped)
RF_Full_S_Avg = mean(tmpRF,4);
RF_Full_S_Avg(:) = RF_Full_S_Avg(:)/abs(min(RF_Full_S_Avg(:)))


tmpRF = [];
for i = 1:numel(lstDTR_S)
    tmpRF(:,:,:,i) = revolveVector3(eval(lstDTR_S{i})', sz, ONrgc);
end
RF_DTR_S_Avg = mean(tmpRF,4);
RF_DTR_S_Avg(:) = RF_DTR_S_Avg(:)/abs(min(RF_DTR_S_Avg(:)));


tmpRF = [];
for i = 1:numel(lstPartial_S)
    tmpRF(:,:,:,i) = revolveVector3(eval(lstPartial_S{i})', sz, ONrgc);
end
RF_Partial_S_Avg = mean(tmpRF,4);
RF_Partial_S_Avg(:) = RF_Partial_S_Avg(:)/abs(min(RF_Partial_S_Avg(:)));


%% Show final sptatio-temporal RF as a movie
tmpRF = RF_Full_S_Avg;
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

%% ON cell analysis
% Selective batch analysis, for cases where you have only one group (sustained or transient)
Kernels = cat(4, RF_Full_S_Avg, RF_Partial_S_Avg, RF_DTR_S_Avg);
groupNames = {'Full_ON_S', 'Partial_ON_S', 'DTR_ON_S'};
[MSEBatch, PSNRBatch, SSIMBatch, diffErr] = analyzeSelectiveBatchGroups(V, Kernels, videofilename, groupNames);

SSIMbatch_avg = mean(SSIMBatch,2)';
MSEbatch_avg  = mean(MSEBatch,2)';  
PSNRbatch_avg = mean(PSNRBatch,2)';
diffErr_avg = mean(diffErr,2)';;
