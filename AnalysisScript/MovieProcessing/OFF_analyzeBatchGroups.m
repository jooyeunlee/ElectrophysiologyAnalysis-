% Written by Luca Della Santina
% Edited by Jeremiah V John

function [MSEbatch, PSNRbatch, SSIMbatch, diffErr, SSIM_filtered] = analyzeBatchGroups(V, Kernels, videofilename, groupNames)
%% Convolve video with each kernel in `Kernels` for each group
%   V (w,h,frames): uint8 encoding of video to convolve
%   Kernels (w,h,temporal, group): Spatiotemporal RF (w,h,t) for each group
%   videofilename (char): filename of video
%   groupNames (1,N): string array for each group (F_S,P_S,D_S,F_T,P_T,D_T)
ksize = size(Kernels);
assert(isequal(ksize(end), length(groupNames)),["Array of Kernels is not equal to number of group names. Kernels length: "  size(Kernels,2) " and groupName length: " size(groupNames,2)])
assert(sum(contains(lower(groupNames), "full"))==2, 'To create image masks between control and other conditions, please include the tags `Full` and `Full` in the group name to specify the two controls RF.')
assert(sum(contains(lower(groupNames), ["_s" "_t"])) >= 2, 'Please label group name with `_S_` or `_T_` for sustained and transient, respectively')
%
sustained = groupNames(contains(lower(groupNames), '_s'));
transient = groupNames(contains(lower(groupNames), '_t'));

% To store filtered videos for each group
for group=1:size(groupNames,2)
    Vfilts.(groupNames{group}) = zeros(size(V,1), size(V,2), size(V,3)-size(Kernels,3));
end
%In this case, metrics are a 2D matrix, where each row is a group
%(Control-S, Partial-S, DTR-S, Control-T, ...)
SSIMbatch = zeros(size(groupNames,2),size(V,3)-size(Kernels,3));
MSEbatch  = zeros(size(groupNames,2),size(V,3)-size(Kernels,3));
PSNRbatch = zeros(size(groupNames,2),size(V,3)-size(Kernels,3));
diffErr = zeros(size(groupNames,2),size(V,3)-size(Kernels,3));

for group = 1:length(groupNames)
    disp(['Starting ' groupNames{group}]);
    [MSEclip, PSNRclip, SSIMclip, Vfilts.(groupNames{group})] = analyzeGroup(V,Kernels(:,:,:,group),videofilename,groupNames{group});
    MSEbatch(group,:) = MSEclip;
    PSNRbatch(group,:) = PSNRclip;
    SSIMbatch(group,:) = SSIMclip;
end
% Create subtraction mask between control and (partial or DTR) for
% sustained
disp('Subtracting control videos with partial and DTR videos')

suFullNm = sustained{contains(lower(sustained),'full')};
suPartialNm = sustained{contains(lower(sustained),'partial')};
suDTRNm = sustained{contains(lower(sustained),'dtr')};
suDiffPartial = imabsdiff(Vfilts.(suFullNm), Vfilts.(suPartialNm));
suDiffDTR = imabsdiff(Vfilts.(suFullNm), Vfilts.(suDTRNm));
suDiffPartial_DTR = imabsdiff(Vfilts.(suPartialNm), Vfilts.(suDTRNm));

% transient
trFullNm = transient{contains(lower(transient),'full')};
trPartialNm = transient{contains(lower(transient),'partial')};
trDTRNm = transient{contains(lower(transient),'dtr')};
trDiffPartial = imabsdiff(Vfilts.(trFullNm), Vfilts.(trPartialNm));
trDiffDTR = imabsdiff(Vfilts.(trFullNm), Vfilts.(trDTRNm));
trDiffPartial_DTR = imabsdiff(Vfilts.(trPartialNm), Vfilts.(trDTRNm));

% % Save masked filtered movies for sustained and transient
% saveDir = 'MoviesProcessed_[71,71]/';
% if ~exist(saveDir, 'dir')
%     mkdir(saveDir)
% end
% fname = [extractBefore(videofilename, '.mp4') '_'];
% saveToVideo(suDiffPartial, [saveDir fname '_S_' '_MaskPartial']);
% saveToVideo(suDiffDTR, [saveDir fname '_S_' '_MaskDTR']);
% saveToVideo(suDiffPartial_DTR, [saveDir fname '_S_' '_MaskPartial_DTR']);
% saveToVideo(trDiffPartial, [saveDir fname '_T_' '_MaskPartial']);
% saveToVideo(trDiffDTR, [saveDir fname '_T_' '_MaskDTR']);
% saveToVideo(trDiffPartial_DTR, [saveDir fname '_T_' '_MaskPartial_DTR']);

% Calculate errors between filtered Control and (Partial or DTR) = MSE
disp('Calculating pixel-wise error between control and (Partial or DTR) videos')
for group=1:length(groupNames) 
    if any(contains(lower(groupNames{group}), '_s'))
        %sustained
        diffErr(group,:) = mean((Vfilts.(groupNames{group}) - Vfilts.(suFullNm)).^2, [1,2]);
    elseif any(contains(lower(groupNames{group}), '_t'))
        %transient
        diffErr(group,:) = mean((Vfilts.(groupNames{group}) - Vfilts.(trFullNm)).^2, [1,2]); 
    end
end
diffErr(length(groupNames)+1,:) = mean((Vfilts.(suPartialNm) - Vfilts.(suDTRNm)).^2, [1,2]); % Adding sustained Partial-DTR (MSE at 7th row)
diffErr(length(groupNames)+2,:) = mean((Vfilts.(trPartialNm) - Vfilts.(trDTRNm)).^2, [1,2]); % Adding transient Partial-DTR (MSE at 8th row)
