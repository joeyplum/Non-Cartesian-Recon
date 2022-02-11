function Image_Out = StandardRecon(AcqSize,data,traj,DCF,weights,ReconSize,PixelShift)
%% A Function written to reconstruct Images when K-space data and trajectories are passed to it
% Uses Pipe's Group DCF and gridding code. This is for 3D data
% 
% AcqSize - Scalar: Output image matrix size
%
% data - KSpace Data in (N_ro x N_Proj) or (N_ro x N_Proj x N_Chan)
%
% traj - point in kspace corresponding to the data vector - columns for
% x,y, and z. (3 x N_ro x N_Proj)
%
% DCF - optional; DCF to use, will skip calculating it
%
% weights - optional; weights for preferentially weighting data, if passed,
%   will calculated new DCF even if one is passed
%
% ReconSize - optional; matrix size to interpolate to
%
% PixelShift - array of size 3 with pixels to shift by

%% Analyze inputs
PreWeight = 0; %not passed
CalculateDCF = 0; %skip calulation
if exist('DCF','var')==0 %if not passed, calculate
    CalculateDCF = 1;
elseif isempty(DCF) %if empty, calculate
    CalculateDCF = 1;
end
if exist('weights','var')==1 %if passed, calculate
    if isempty(weights)
        %if empty, skip
    else
        PreWeight = 1;
    end
end
if exist('ReconSize','var')==0 %if not passed, set to no interp
    ReconSize = AcqSize;
elseif isempty(ReconSize) %if empty, set to no interp
    ReconSize = AcqSize;
end
if exist('PixelShift','var')==0 %if not passed, set to 0's
    PixelShift = [0, 0, 0];
elseif isempty(PixelShift) %if empty, set to 0's
    PixelShift = [0, 0, 0];
end

%% Settings
numIter = 5;
effMtx  = AcqSize;
alpha = 2; %Image overgrid factor
osf = 2.1; %grid oversample factor on top of effMtx, at least 2, 2.1 optimal, >2.6 conservative
verbose = 0; %c verbose is not real time so not very useful
numThread = 1; %unthreaded since not Posix

%% Prep inputs
data_real = reshape(real(data), [1 size(data)]);
data_imag = reshape(imag(data), [1 size(data)]);
data_comb = cat(1, data_real, data_imag);
data_comb = double(data_comb);
traj = double(traj);
% weights = double(weights);

clear data data_real data_imag

%% DCF
if (PreWeight)
    disp('Calculating weighted DCF...');
    DCF = recon.DC.sdc3_MAT(traj,numIter,effMtx,verbose,osf,weights);
    DCF = double(DCF);
    disp('Calculated weighted DCF.');
elseif (CalculateDCF)
    disp('Calculating DCF...');
    DCF = recon.DC.sdc3_MAT(traj,numIter,effMtx,verbose,osf);
    DCF = double(DCF);
    disp('Calculated DCF.');
else
    DCF = double(DCF);
end

%% Grid
gdata = zeros(2,effMtx*alpha,effMtx*alpha,effMtx*alpha,size(data_comb,4));
parfor coil = 1:size(data_comb,4)
    disp(['Gridding channel ', num2str(coil), ' of ', num2str(size(data_comb,4)), '...']);
    gdata(:,:,:,:,coil) = recon.Grid.grid3_MAT(squeeze(data_comb(:,:,:,coil)),traj,DCF,effMtx*alpha,numThread);
end
gdata = squeeze(gdata(1,:,:,:,:) + 1j*gdata(2,:,:,:,:));

clear data_comb traj DCF

%% FFT to image
disp('FFTing gridded data...');
Image = fft(gdata,ReconSize*alpha,1);
Image = fft(Image,ReconSize*alpha,2);
Image = fft(Image,ReconSize*alpha,3);
Image = fftshift(Image,1);
Image = fftshift(Image,2);
Image = fftshift(Image,3);

clear gdata

%% Calculated rolloff kernel
disp('Calculating rolloff kernel...');
delta = [1.0, 0.0];
k_not = [0.0, 0.0, 0.0];
DCF_not = 1.0;
rokern = recon.Grid.grid3_MAT(delta',k_not',DCF_not,effMtx*alpha,numThread);

%% FFT to rolloff image
disp('FFTing rolloff kernel to image space...');
rokern = squeeze(rokern(1,:,:,:) + 1j*rokern(2,:,:,:));
rokern = fftn(rokern, [ReconSize*alpha ReconSize*alpha ReconSize*alpha]);
rokern = fftshift(rokern,1);
rokern = fftshift(rokern,2);
rokern = fftshift(rokern,3);

%% Apply rolloff kernel image
for coil = 1:size(Image,4)
    disp(['Applying rolloff for channel ', num2str(coil), ' of ', num2str(size(Image,4)), '...']);
    Image(:,:,:,coil) = Image(:,:,:,coil) ./ rokern;
end

%% Shift pixels
disp('Applying shift and cropping image...');
Image = circshift(circshift(circshift(Image,round(PixelShift(1)),1),round(PixelShift(2)),2),round(PixelShift(3)),3);

%% Crop
xs = floor(ReconSize - ReconSize/alpha)+1;
xe = floor(ReconSize + ReconSize/alpha);
Image_Out = Image(xs:xe,xs:xe,xs:xe,:);

disp('Reconstruction complete.');

end





