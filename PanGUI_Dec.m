%% Image Deconvolution:
% Deconvolving the video (unblurring the Point Spread Function - PSF) into the Tissue Reflectivity Function (TRF).
%
% Syntax:
% shortvid_ceus = PanGUI_Dec(shortvid,name,fs)
%
% Input:
% shortvid - Original second video (CEUS) double[0,255].
% name - The file's name.
% fs - The video's frame rate.
%
% Output:
% shortvid_ceus - Registration of the second video (CEUS) double[0,255].

function shortvid_ceus = PanGUI_Dec(shortvid,name,fs)
% Results directory:
respar = '.\Results - temp\';

% Changing the number of frames in the analysis:
numofframes = size(shortvid,4); % Default: defined in the parameters.

% Choosing an image:
framestart = 1;
IMG0 = im2double(uint8(shortvid(:,:,:,framestart)));
%IMG0 = shortvid(:,:,:,framestart);

% Choose M0 and N0 automatically:
L_m = size(IMG0,1);
L_n = size(IMG0,2);
pks = zeros(L_m,L_n);
locs = zeros(L_m,L_n);
m_temp = zeros(L_m,1);
for i = 1:L_m
    [pks_temp,locs_temp] = findpeaks(IMG0(i,:));
    pks(i,1:length(pks_temp)) = pks_temp;
    locs(i,1:length(locs_temp)) = locs_temp;
    test = locs(i,:).*(pks(i,:)>mean(IMG0,'all'));
    test(test == 0) = L_m;
    m_temp(i) = min(test);
end
[N0_d,M0_d] = min(m_temp);
clear pks pks_temp locs locs_temp test;
if N0_d > 5, N0_d = N0_d-5; else, N0_d = 1; end
if M0_d > 5, M0_d = M0_d-5; else, M0_d = 1; end
%
%M0_d = 55; % 1st Row
%N0_d = 85; % 1st Column
M = 10; % # of Rows
N = 10; % # of Columns
%IMG1 = cat(3,IMG0,IMG0,IMG0);
IMG1 = IMG0;
IMG1(M0_d:M0_d+M-1,N0_d,1) = 1;
IMG1(M0_d:M0_d+M-1,N0_d+N-1,1) = 1;
IMG1(M0_d,N0_d:N0_d+N-1,1) = 1;
IMG1(M0_d+M-1,N0_d:N0_d+N-1,1) = 1;

% Choosing an image part:
IMG = shortvid(M0_d:M0_d+M-1,N0_d:N0_d+N-1,1,framestart);

% Plotting the frame:
if ~exist([respar,'Deconvolution\'],'dir')
    mkdir([respar,'Deconvolution\']);
end
% Original Framed:
f = figure('visible','off');
imagesc(IMG0);
colormap gray;
title('Original');
saveas(gcf,[respar,'Deconvolution\Original.jpg']);
close(f);
% Original Framed:
f = figure('visible','off');
imagesc(IMG1);
colormap gray;
title('Original (Kernel in red)');
saveas(gcf,[respar,'Deconvolution\OriginalFramed.jpg']);
close(f);
% Original Focused:
f = figure('visible','off');
imagesc(IMG);
colormap gray;
title('Kernel');
saveas(gcf,[respar,'Deconvolution\Kernel.jpg']);
close(f);

% Looking for Kernel - version 2 (by deconvolution properties):
max_dif = 1e-10;
max_loop = 20;
change_delay = 0;
type = 1;
nn = 4;
param_struct0 = {max_dif,max_loop,change_delay,type,nn};
[kerfit,sigma_m,sigma_n] = ker_check_2(IMG,param_struct0); % IMG - short, IMG0 - full.
% Only for the Phantom:
% sigma_m = 0.7;
% sigma_n = 0.7;

% The de-convolution algorithm:
max_loop = 60;
change_delay = fix(max_loop/15); % "fix(iter/15)" -> accuracy of ~0.01. Was 3.
win = 0; % '1' - Include Hann window, '0' - without window, '-1' for Sphere (ellipsoid) Kernel.
param_struct = {max_dif,max_loop,change_delay,type,nn};
if win == -1 % Ellipsoid
    [Xmax,Xavg,Yrf_mat,Yrec_mat,Yrfmax,Yrfavg] = decon5_p3(vid,numofframes,2*sigma_m,2*sigma_n,param_struct,win);
else
    [Xmax,Xavg,Yrf_mat,Yrec_mat,Yrfmax,Yrfavg] = decon5_p3(shortvid,numofframes,sigma_m,sigma_n,param_struct,win);
end
Yrfmax_thr = 0.05;
gamma = 0.7;

% Plot All Frames:
% The de-convolution algorithm version to use:
% Default: 'RegCons', Seems to bring the best result.
ver = 'RegCons';

% Total RF sum scaled
f = figure('visible','off');
imagesc(Yrfmax);
colormap gray;
title(['Yrfmax ''RegCons'', \sigma_m = ',num2str(sigma_m),', \sigma_n = ',num2str(sigma_n),', scaled, ',int2str(numofframes),' frames']);
saveas(gcf,[respar,'Deconvolution\Yrfmax_',ver,'_sc.jpg']);
close(f);
%
f = figure('visible','off');
imagesc(Yrfmax.*(Yrfmax > Yrfmax_thr));
colormap gray;
title(['Yrfmax ''RegCons'', \sigma_m = ',num2str(sigma_m),', \sigma_n = ',num2str(sigma_n),', thr = ',num2str(Yrfmax_thr),', scaled, ',int2str(numofframes),' frames']);
saveas(gcf,[respar,'Deconvolution\Yrfmax_thr_',ver,'_sc.jpg']);
close(f);

% Total RF average scaled
f = figure('visible','off');
imagesc(Yrfavg);
colormap gray;
title(['Yrfavg ''RegCons'', \sigma_m = ',num2str(sigma_m),', \sigma_n = ',num2str(sigma_n),', scaled, ',int2str(numofframes),' frames']);
saveas(gcf,[respar,'Deconvolution\Yrfavg_',ver,'_sc.jpg']);
close(f);

% Total RF average scaled with gamma correction
f = figure('visible','off');
imagesc((Yrfavg/max(Yrfavg(:))).^gamma);
colormap gray;
title(['Yrfavg ''RegCons'', \sigma_m = ',num2str(sigma_m),', \sigma_n = ',num2str(sigma_n),', \gamma = ',num2str(gamma),', scaled, ',int2str(numofframes),' frames']);
saveas(gcf,[respar,'Deconvolution\Yrfavg_',ver,'_sc.jpg']);%
close(f);

% Maximum Intensity Persistence scaled
f = figure('visible','off');
imagesc(Xmax);
colormap gray;
title(['Xmax - MIP, scaled, ',int2str(numofframes),' frames']);
saveas(gcf,[respar,'Deconvolution\Xmax_sc.jpg']);
close(f);

% Total RF average scaled
f = figure('visible','off');
imagesc(Xavg);
colormap gray;
title(['Xavg - Mean, scaled, ',int2str(numofframes),' frames']);
saveas(gcf,[respar,'Deconvolution\Xavg_sc.jpg']);
close(f);

% Reconstruction:
f = figure('visible','off');
imagesc(Yrec_mat(:,:,1));
colormap gray;
title('Reconstruction - frame 1');
saveas(gcf,[respar,'Deconvolution\Reconstruction.jpg']);
close(f);

% Reflectivity function:
f = figure('visible','off');
imagesc(Yrf_mat(:,:,1));
colormap gray;
title('Reflectivity function - frame 1');
saveas(gcf,[respar,'Deconvolution\Reflectivity.jpg']);
close(f);

% Saves the data:
save([respar,'Deconvolution\decon.mat'],'param_struct0','param_struct','Xmax','Xavg','Yrf_mat','Yrec_mat','Yrfmax','Yrfavg');


% Preparing for "VideoWriter".
% Deconvolved:
%Yrf_mat_out = Yrf_mat;
Yrf_mat_out = zeros(size(Yrf_mat,1),size(Yrf_mat,2),3,size(Yrf_mat,3));
%Yrec_mat_out = Yrec_mat;
Yrec_mat_out = zeros(size(Yrec_mat,1),size(Yrec_mat,2),3,size(Yrec_mat,3));
for i = 1:numofframes
    Yrf_mat_out(:,:,1,i) = 255*Yrf_mat(:,:,i);
    Yrf_mat_out(:,:,2,i) = 255*Yrf_mat(:,:,i);
    Yrf_mat_out(:,:,3,i) = 255*Yrf_mat(:,:,i);
    Yrec_mat_out(:,:,1,i) = 255*Yrec_mat(:,:,i);
    Yrec_mat_out(:,:,2,i) = 255*Yrec_mat(:,:,i);
    Yrec_mat_out(:,:,3,i) = 255*Yrec_mat(:,:,i);
end

% Saving the uncompressed video:
shortvid_ceus = Yrf_mat_out;

%{
Yrf_mat_out = uint8(Yrf_mat_out);
% Intensity correction:
%Yrf_mat_out = uint8(cat(3,Yrf_mat_out,cat(3,Yrf_mat_out,Yrf_mat_out)));

% Saves the video:
% Deconvolved:
%v = VideoWriter([respar,'Deconvolved - ',name],'Grayscale AVI'); % Saturated video
v = VideoWriter([respar,'Deconvolved_CEUS - ',name],'Uncompressed AVI'); % Saturated video
v.FrameRate = fs;
open(v);
writeVideo(v,Yrf_mat_out);
close(v);
%}

% Saves Original videos:
save([respar,'Deconvolved_CEUS - ',name,'.mat'],'shortvid_ceus');
save([respar,'Reconstructed_CEUS - ',name,'.mat'],'Yrec_mat_out');

% Saving memory space:
clear v;
clear Yrf_mat_out Yrec_mat_out;
clear Yrec_mat;
end