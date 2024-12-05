% Written by Luca Della Santina
% Edited by Jeremiah V John

function [MSEclip, PSNRclip, SSIMclip, filteredV] = analyzeGroup(V, Kernel, videofilename, groupName)
%% Convolve the array V with Kernel, then analyze.
disp('Convolving movie with filter using FFT...');
tic
Fast_VP = dftConv3(V, Kernel);

disp(['Done in...' num2str(toc) ' seconds']);

VP = Fast_VP;
SSIMclip = zeros(1,size(V,3)-size(Kernel,3));
MSEclip  = zeros(1,size(V,3)-size(Kernel,3));
PSNRclip = zeros(1,size(V,3)-size(Kernel,3));

% 202030825 to create 3D array for original and filtered video(JJ)

originalV = zeros(size(V,1), size(V,2), size(V,3)-size(Kernel,3));
filteredV = zeros(size(V,1), size(V,2), size(V,3)-size(Kernel,3));
% diffV = zeros(size(V,1), size(V,2), size(V,3)-size(Kernel,3));

for frame = 1:size(V,3)-size(Kernel,3)
    % Convert frame to 8bit (should not be necessary if input V is 8bit)
    frameV = V(:,:,frame);
    frameV = frameV - min(frameV);
    frameV8bit = uint8(round(255*(frameV/max(frameV(:)))));
    originalV(:,:,frame) = frameV;              % 20230825
    
    % Convert each processed frame to 8bit (since convn outputs double)
    frameVP = VP(:,:,frame);
    frameVP = frameVP - min(frameVP(:));
    frameVP = uint8(round(255*(frameVP/max(frameVP(:)))));
    filteredV(:,:,frame) = frameVP;             % 20230825
    
%     %    imshowpair(frameV, frameVP, 'montage');
    MSEclip(frame)  = immse(frameVP, frameV);
    PSNRclip(frame) = psnr(frameVP, frameV);
    SSIMclip(frame) = ssim(frameVP, frameV);
% %     diff_clip = imabsdiff(frameV, frameVP);          % can use abs(x-y)
% %     diff_frame(frame) = mean((frameVP - frameVP_dtr).^2, 'all'); % mean square the difference
% %     diffV(:,:,frame) = diff_clip;
end

%Save greyscale original and filtered movies.
saveDir = 'MoviesProcessed_ON_[83,83]/'; % change folder name based on cell types and size
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end
fname = [extractBefore(videofilename, '.mp4') '_' groupName];
saveToVideo(originalV, [saveDir fname '_original']);
saveToVideo(filteredV, [saveDir fname '_filtered']);
% saveToVideo(diffV, ['MoviesProcessed/' fname '_diff']);

