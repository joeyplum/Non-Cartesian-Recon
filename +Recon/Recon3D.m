function [Image_Out, DCF] = Recon3D(AcqSize,data,traj,DCF,weights,ReconSize,PixelShift)
%% A Function written to reconstruct Images when K-space data and trajectories are passed to it
% This code uses Pipe's Group DCF and gridding code. 
% The code will take 3D k-space data, and return a 3D image. 
%
%
% Inputs:
% 
% AcqSize - single scalar value = output image matrix size. 
% Example: if you acquire an image with a 128*128 matrix, enter 128 and you
% will reconstruct a full image that is 128 pixels wide. If you enter 200,
% for example, then the reconstructed image will be 200 pixels wide, and
% you will reconstruct data outside of your aquired FOV. 
%
% data - KSpace Data in (N_samples x N_excitations) or (N_samples x
% N_excitations x N_Chan).
%
% traj - Coordinate locations in kspace corresponding to the data vector:
% columns for x,y, and z (3 x N_samples x N_excitations). 
%
% Additional (optional) parameters:
% DCF - Use a home-made DCF instead of calculating one. 
%
% weights - Weights for preferentially weighting data. If passed, a new DCF
% will be calculated, even if a DCF has been passed.
%
% ReconSize - Matrix size to interpolate to. This feature is pretty neat if
% you want to for example: downsample a 200*200 image to 100*100 size, or
% vice versa.
%
% PixelShift - Array of size 3 to shift pixels by ([shift_x, shift_y, shift_z]).
%
%
% Outputs:
%
% Image_Out - 2D image (complex) of the reconstructed k-space data.
%
% DCF - export the DCF for your k-space trajectories. The exported DCF will
% be of size (N_samples x N_excitations).
%
%
% Example usage for command line:
% figure; [Image_Out, DCF] = Recon.Recon3D(116,data,traj,[],[],116,[]); imagesc(abs(Image_Out(:,:,1)))

%% Initialize and analyze inputs:
CalculateDCF = 0; % Initialize to skip calculation.
if exist('DCF','var') == 0 % If no DCF is passed, change to calculate DCF.
    CalculateDCF = 1;
elseif isempty(DCF) % If a DCF is set to [] (empty), then calculate DCF.
    CalculateDCF = 1;
end

PreWeight = 0; % Initialize to not pass.
if exist('weights','var') == 1 % If weights are passed, then calculate.
    if isempty(weights)
        % If empty, then skip.
    else
        PreWeight = 1;
        % Set PreWeight to calculate.
    end
else 
    weights = [];
end

if exist('ReconSize','var') == 0 % If not passed, do not interpolate image.
    ReconSize = AcqSize;
elseif isempty(ReconSize) % If empty, do not interpolate image.
    ReconSize = AcqSize;
end

if exist('PixelShift','var') == 0 % If not passed, set to zero pixel shift.
    PixelShift = [0, 0, 0];
elseif isempty(PixelShift) % If empty, set to zero pixel shift.
    PixelShift = [0, 0, 0];
end

%% Reconstruction settings:
numIter = 5;
effMtx  = AcqSize;
alpha = 2; % Image overgrid factor.
osf = 2.1; % Grid oversample factor on top of effMtx.
% Recommended: at least 2, 2.1 optimal, >2.6 conservative.
verbose = 0; 
% Recommended: 0, as c verbose is not real time so not very useful.
numThread = 1; % Unthreaded since not Posix.

%% Prepare inputs:
data_real = reshape(real(data), [1 size(data)]);
data_imag = reshape(imag(data), [1 size(data)]);
data_comb = cat(1, data_real, data_imag);
data_comb = double(data_comb);
traj = double(traj);
weights = double(weights);

clear data data_real data_imag

%% Calculate DCF:
if (PreWeight)
    disp('Calculating weighted DCF...');
    DCF = Recon.DC.sdc3_MAT(traj,numIter,effMtx,verbose,osf,weights);
    DCF = double(DCF);
    disp('Calculated weighted DCF.');
elseif (CalculateDCF)
    disp('Calculating DCF...');
    DCF = Recon.DC.sdc3_MAT(traj,numIter,effMtx,verbose,osf);
    DCF = double(DCF);
    disp('Calculated DCF.');
else
    DCF = double(DCF);
end

%% Grid data:
gdata = zeros(2,effMtx*alpha,effMtx*alpha,effMtx*alpha,size(data_comb,4));
parfor coil = 1:size(data_comb,4)
    disp(['Gridding channel ', num2str(coil), ' of ', num2str(size(data_comb,4)), '...']);
    gdata(:,:,:,:,1,coil) = Recon.Grid.grid3_MAT(squeeze(data_comb(:,:,:,coil)),traj,DCF,effMtx*alpha,numThread);
end
gdata = squeeze(gdata(1,:,:,:,:) + 1j*gdata(2,:,:,:,:));
disp('Completed gridding.')

clear data_comb traj 

%% FFT to image:
disp('FFTing gridded data...');
Image = fft(gdata,ReconSize*alpha,1);
Image = fft(Image,ReconSize*alpha,2);
Image = fft(Image,ReconSize*alpha,3);
Image = fftshift(Image,1);
Image = fftshift(Image,2);
Image = fftshift(Image,3);
disp('Completed FFT of gridded data.');

clear gdata

%% Calculated rolloff kernel:
disp('Calculating rolloff kernel...');
delta = [1.0, 0.0];
k_not = [0.0, 0.0, 0.0];
DCF_not = 1.0;
rokern = Recon.Grid.grid3_MAT(delta',k_not',DCF_not,effMtx*alpha,numThread);
disp('Completed calculation of rolloff kernel.');

%% FFT to rolloff image:
disp('FFTing rolloff kernel to image space...');
rokern = squeeze(rokern(1,:,:,:) + 1j*rokern(2,:,:,:));
rokern = fftn(rokern, [ReconSize*alpha ReconSize*alpha ReconSize*alpha]);
rokern = fftshift(rokern,1);
rokern = fftshift(rokern,2);
rokern = fftshift(rokern,3);
disp('Completed FFT of rolloff kernel to image space.');

%% Apply rolloff kernel image:
for coil = 1:size(Image,4)
    disp(['Applying rolloff for channel ', num2str(coil), ' of ', num2str(size(Image,4)), '...']);
    Image(:,:,:,coil) = Image(:,:,:,coil) ./ rokern;
end

%% Shift pixels:
disp('Applying shift and cropping image...');
Image = circshift(circshift(circshift(Image,round(PixelShift(1)),1),round(PixelShift(2)),2),round(PixelShift(3)),3);

%% Crop:
xs = floor(ReconSize - ReconSize/alpha)+1;
xe = floor(ReconSize + ReconSize/alpha);
Image_Out = Image(xs:xe,xs:xe,xs:xe,:);

disp('Reconstruction complete.');

end


