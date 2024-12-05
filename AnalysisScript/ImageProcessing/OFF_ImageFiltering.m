% Written by Luca Della Santina

%% Accumulate spatio-temporal RF per condition
lstFull = {'space0', 'space1', 'space2', 'space3', 'space4', 'space5', 'space6', 'space7', 'space8', 'space9', 'space10', 'space11', 'space12'};%, 'space13', 'space14', 'space15', 'space16'};
for i = 1:numel(lstFull)
    [RFfull(:,:,i), RFfullSigned(:,:,i)] = revolveVector_OFF(eval(lstFull{i}));
end
RFfullAvg = mean(RFfull,3);
RFfullSignedAvg = mean(RFfullSigned,3);

lstHalf = {'subspace0', 'subspace1', 'subspace2', 'subspace3', 'subspace4', 'subspace5', 'subspace6', 'subspace7', 'subspace8', 'subspace9', 'subspace10', 'subspace11', 'subspace12'};%,  'subspace13', 'subspace14', 'subspace15',  'subspace16'};
for i = 1:numel(lstHalf)
    [RFhalf(:,:,i), RFhalfSigned(:,:,i)] = revolveVector_OFF(eval(lstHalf{i}));
end
RFhalfAvg = mean(RFhalf,3);
RFhalfSignedAvg = mean(RFhalfSigned,3);


lstDTR = {'dtrspace0', 'dtrspace1', 'dtrspace2', 'dtrspace3', 'dtrspace4', 'dtrspace5', 'dtrspace6', 'dtrspace7', 'dtrspace8', 'dtrspace9', 'dtrspace10', 'dtrspace11', 'dtrspace12', 'dtrspace13', 'dtrspace14', 'dtrspace15'};
lstDTRflip = {false, true, false, true, false, false, false, false, false, false, false, false}; % Need to match polarity for each cell (false for blue center, true for red center to flip)
for i = 1:numel(lstDTR)
    [RFdtr(:,:,i), RFdtrSigned(:,:,i)]  = revolveVector_OFF(eval(lstDTR{i}));
end
RFdtrAvg = mean(RFdtr,3);
RFdtrSignedAvg = mean(RFdtrSigned,3);

%% Plot Average spatio-temporal RF and calculate SSIM
figure('Units', 'normalized', 'Position', [0.1, 0, 0.8, 0.9]);
tiledlayout(2,3, 'Tilespacing' , 'tight');
nexttile;
%imshow(uint8(round(255*RFfullAvg)));
imagesc(RFfullSignedAvg);
colormap(jet);
title('Full-stim');
axis off;

nexttile;
%imshow(uint8(round(255*RFhalfAvg)));
imagesc(RFhalfSignedAvg);
colormap(jet);
title('Half-stim');
axis off;

nexttile;
%imshow(uint8(round(255*RFdtrAvg)));
imagesc(RFdtrSignedAvg);
colormap(jet);
colorbar;
title('DTR');
axis off;

nexttile;
%plot(RFfullAvg(ceil(size(RFfullAvg,1)/2),:), '-k');
%hold on;
spacescale = -1:1/floor(size(RFfullSignedAvg,1)/2):1;
plot(spacescale, RFfullSignedAvg(ceil(size(RFfullSignedAvg,1)/2),:), '-r');
ylabel('Magnitude (a.u.)');

nexttile;
%plot(RFhalfAvg(ceil(size(RFhalfAvg,1)/2),:), '-k');
%hold on;
spacescale = -1:1/floor(size(RFhalfSignedAvg,1)/2):1;
plot(spacescale, RFhalfSignedAvg(ceil(size(RFhalfSignedAvg,1)/2),:), '-r');
xlabel('Space (mm)');

nexttile;
%plot(RFdtrAvg(ceil(size(RFdtrAvg,1)/2),:), '-k');
%hold on;
spacescale = -1:1/floor(size(RFdtrSignedAvg,1)/2):1;
plot(spacescale, RFdtrSignedAvg(ceil(size(RFdtrSignedAvg,1)/2),:), '-r');

clear i lst* RF* ans I* Similarity fName