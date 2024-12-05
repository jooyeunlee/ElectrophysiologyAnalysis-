% Before running this script, make sure to have the variables 
% RFfullSignedAvg, RFhalfSignedAvg, RFdtrSignedAvg in matlab workspace
%
% Version 1.1 
%   Stored max of SSIM obtained from and inverse of the photo
%   Added 0.5 as DownscaleFactor
% Version 1.2
%   Verbose progress of current image number
%   Added calculation of PSNR and RMSE
% Version 1.3
%   Added calculation of SNR
% Version 1.4
%   Manual calculation of PSNR from MSE

fPath = uigetdir;
files = dir(fPath);
files = files(~[files.isdir]);

DownscaleFactors = [0.5, 1, 2, 5, 10, 50, 100, 500];
%DownscaleFactors = [2, 5, 10, 50, 100, 500];
IrfRatios = DownscaleFactors*1/2000;
disp(['Image/RF ratios: ' num2str(IrfRatios)]);
disp(' ');

% Structural Similarity Index
SSIMfull = zeros(numel(DownscaleFactors),numel(files));
SSIMhalf = zeros(numel(DownscaleFactors),numel(files));
SSIMdtr = zeros(numel(DownscaleFactors),numel(files));

% Peak signal-to-noise ratio
PSNRfull = zeros(numel(DownscaleFactors),numel(files));
PSNRhalf = zeros(numel(DownscaleFactors),numel(files));
PSNRdtr = zeros(numel(DownscaleFactors),numel(files));

% Mean square error
MSEfull = zeros(numel(DownscaleFactors),numel(files));
MSEhalf = zeros(numel(DownscaleFactors),numel(files));
MSEdtr = zeros(numel(DownscaleFactors),numel(files));


for i = 1:numel(files)
    % Convolve Natural Images and create Image/Filter vs SSIM plot
    
    fName = [files(i).folder filesep files(i).name];    
    disp(['Loading photo ' num2str(i) '/' num2str(numel(files)) ': ' fName]);
    I = imread(fName);
    Ibw = rgb2gray(I);
    
        
    for ds = 1:numel(DownscaleFactors)
        disp(['  |- Downsampling by ' num2str(DownscaleFactors(ds))]);
        IbwDS = imresize(Ibw, 1/DownscaleFactors(ds), 'nearest');
        
        Ifull = conv2(IbwDS, RFfullSignedAvg,'same');
        Ifull8bit = Ifull - min(Ifull(:));
        Ifull8bit = uint8(255*(round(Ifull8bit/max(Ifull8bit(:)))));
        SSIMfull(ds, i) = max([ssim(Ifull8bit, IbwDS), ssim(255-Ifull8bit, IbwDS)]);
        MSEfull(ds, i) = max([immse(Ifull8bit, IbwDS), immse(255-Ifull8bit, IbwDS)]);
        PSNRfull(ds, i) = 10*log10(255^2/MSEfull(ds,i));
        
        Ihalf = conv2(IbwDS,RFhalfSignedAvg,'same');
        Ihalf8bit = Ihalf - min(Ihalf(:));
        Ihalf8bit = uint8(255*(round(Ihalf8bit/max(Ihalf8bit(:)))));
        SSIMhalf(ds, i) = max([ssim(Ihalf8bit, IbwDS), ssim(255-Ihalf8bit, IbwDS)]);
        MSEhalf(ds, i) = max([immse(Ihalf8bit, IbwDS), immse(255-Ihalf8bit, IbwDS)]);
        PSNRhalf(ds, i) = 10*log10(255^2/MSEhalf(ds,i));
        
        Idtr = conv2(IbwDS,RFdtrSignedAvg,'same');
        Idtr8bit = Idtr - min(Idtr(:));
        Idtr8bit = uint8(255*(round(Idtr8bit/max(Idtr8bit(:)))));
        SSIMdtr(ds, i) = max([ssim(Idtr8bit, IbwDS),ssim(255-Idtr8bit, IbwDS)]);        
        MSEdtr(ds, i) = max([immse(Idtr8bit, IbwDS),immse(255-Idtr8bit, IbwDS)]);  
        PSNRdtr(ds, i) = 10*log10(255^2/MSEdtr(ds,i));
        
    end        
end