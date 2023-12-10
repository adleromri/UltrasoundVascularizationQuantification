%% De-Convolution Function
% Syntax:
% [Xmax,Xavg,Yrf_mat,Yrec_mat,Yrfmax,Yrfavg] = decon5_p2(vid,numofframes,m_sig2,n_sig2,param_struct,win)
%
% Input:
% vid - the video to analyze. (Double or VideoReader)
% numofframes (optional) - number of frames to analyze.
% m_sig2 - the axial (y-axis) sigma (squared) of the gaussian kernel.
% n_sig2 - the lateral (x-axis) sigma (squared) of the gaussian kernel.
% param_struct - parameters structure, including: max_dif, max_loop,
% change_delay, type, nn.
% win - '1' includes Hann window, '0' without any window, '-1' for Sphere (ellipsoid) Kernel.
%
% Other parameters:
% state (optional) - initial state: 'thr' - soft threshold, 'zero'.
% ver (optional) - version:
% 'HalfCons','DecCons','SevCons','MaxCons','MoreCons','RegCons','DifVal'.
% fil (optional) - filtering: 'on','off'.
% sat (optional) - saturaion: 'on','off'.
%
% Output:
% Xmax - MIP of the original frames.
% Xavg - mean value of the original frames.
% Yrf_mat - the deconvolved images.
% Yrec_mat - the reconstructed images.
% Yrfmax - MIP of the deconvolved images.
% Yrfavg - mean value of the deconvolved images.

function [Xmax,Xavg,Yrf_mat,Yrec_mat,Yrfmax,Yrfavg] = decon5_p3(vid,numofframes,m_sig2,n_sig2,param_struct,win)
%% Loading the data:
if class(vid) == "double"
    if size(vid,4) > 0
        frame0 = vid(:,:,:,1);
        if (size(frame0,3) > 1)&&(max(frame0(:)) > 1)
            X1gray0 = im2double(rgb2gray(uint8(frame0))); % double [0,1]
        elseif size(frame0,3) > 1
            X1gray0 = im2double(rgb2gray(frame0)); % double [0,1]
        elseif max(frame0(:)) > 1
            X1gray0 = im2double(uint8(frame0)); % double [0,1]
        else
            X1gray0 = frame0;
        end
    elseif size(vid,3) > 0
        frame0 = vid(:,:,1);
        if max(frame0(:)) > 1
            X1gray0 = im2double(uint8(frame0)); % double [0,1]
        else
            X1gray0 = frame0; % double [0,1]
        end
    end
elseif class(vid) == "VideoReader"
    frame0 = read(vid,1);
    if size(frame0,3) == 1
        X1gray0 = im2double(frame0); % double [0,1]
    else
        X1gray0 = im2double(rgb2gray(frame0)); % double [0,1]
    end
else
    error('Wrong class input.\n');
end

%% The de-convolution algorithm version to use:
% Need to set "state" (initial state: 'thr' - soft threshold, 'zero').
if ~exist('state','var')
    state = 'zero'; % Default 'zero' (converges a bit more accurate), but 'thr' converges faster.
end
% Need to set "ver" (version: 'HalfCons','DecCons','SevCons','MaxCons','MoreCons','RegCons','DifVal').
%{
if ~exist('ver','var')
    ver = 'RegCons'; % Default. Seems to bring the best result.
end
%}
% Need to set "fil" (filtering: 'on','off').
if ~exist('fil','var')
    fil = 'off'; % Default 'on' (filters noises), but 'off' calculates faster.
end
% Need to set "sat" (saturaion: 'on','off').
if ~exist('sat','var')
    sat = 'off'; % Default.
end

%% Reshaping parameters:
Rmax = size(X1gray0,1);

%% PSF Kernels creation:
% Ammount of PSFs:
N_PSF = 1; % Default: 11 for sectorial, 3 for linear.
% PSFs positions:
y_PSF = Rmax/N_PSF;

%% Gaussian fit:
% Gausian: 1/sqrt(2*pi*sigma^2)*exp(-0.5*(x-mu)^2/sigma^2);
% Log-Gausian: -0.5*(x-mu)^2/sigma^2-0.5*log(2*pi*sigma^2) = -0.5*((x-mu)^2/sigma^2+2*log(sigma))-0.5*log(2*pi);
%max_gauss = 1; % Default: max(calImg(:));
if ~exist('m_sig2','var')
    m_sig2 = 0.8; % 0.8, 0.3, 0.5, Default: 1.5
end
if ~exist('n_sig2','var')
    n_sig2 = 0.8; % 0.8, 0.3, 0.5, Default: 1.5
end
max_gauss = 1/sqrt(2*pi*m_sig2*n_sig2); % Gaussian maximum.
im_g = zeros(size(X1gray0));
mu_g = zeros(2,N_PSF);
for i_mu = 1:N_PSF
    mu_g(:,i_mu) = round([(i_mu-0.5)*y_PSF;0.5*size(im_g,2)]);
end
%
for m = 1:size(im_g,1)
    for n = 1: size(im_g,2)
        coord_g = [m;n];
        % "Gaussian" with top value at 1:
        for i_mu = 1:N_PSF
            % 1 common 2D Gaussian:
            sigma_a = m_sig2;
            sigma_l = n_sig2;
            sig_g = [sigma_a 0;0 sigma_l];
            im_g(m,n) = im_g(m,n) + max_gauss*exp(-0.5*(coord_g-mu_g(:,i_mu)).'*((sig_g)\(coord_g-mu_g(:,i_mu))));
        end
        
        % Normalized gaussian:
        %im_g(m,n) = 1/sqrt((2*pi)^2*det(sig_g)) * exp(-0.5*(coord_g-mu_g).'*((sig_g)\(coord_g-mu_g)));
    end
end
im_g = 1/sum(im_g(:))*im_g; % Normalized to sum 1.
% save('shapeTrans\shapeTrans3\matrices.mat','im_g');
% load('shapeTrans\shapeTrans3\matrices.mat','im_g');

% Sphere (ellipsoid):
if win == -1
    for m2 = 1:size(im_g,1)
        for n2 = 1: size(im_g,2)
            % "Sphere" (ellipsoid) with top value at 1:
            for i_mu = 1:N_PSF
                % top of ellipsoid:
                sigma_a = m_sig2;
                sigma_l = n_sig2;
                im_g(m2,n2) = im_g(m2,n2) + max_gauss*...
                    max(1 - 1/sigma_a*abs(m2 - mu_g(1,i_mu)),0) * max(1 - 1/sigma_l*abs(n2 - mu_g(2,i_mu)),0);
            end
        end
    end
end

%% No reshaping:
Y_im_g = im_g;

%% PSF matrix improvement:
psf_side_m = 7; % Default: 7 or 5.
psf_side_n = 7; % Default: 7 or 5.
mini_side_m = (psf_side_m-1)/2;
mini_side_n = (psf_side_n-1)/2;
psf_mat = zeros(N_PSF*psf_side_m,psf_side_n);
psf_size = [psf_side_m,psf_side_n];
if win
    han = zeros(psf_side_m,psf_side_n);
    for i_m = 1:psf_side_m
        for i_n = 1:psf_side_n
            han(i_m,i_n) = (sin(pi*(i_m-1)/(psf_side_m-1)))^2 * (sin(pi*(i_n-1)/(psf_side_n-1)))^2;
        end
    end
end
for i_psf = 1:N_PSF
    psf_mat(((i_psf-1)*psf_side_m+1):(i_psf*psf_side_m),1:psf_side_n) = ...
        Y_im_g((round((i_psf-0.5)*y_PSF)-mini_side_m):(round((i_psf-0.5)*y_PSF)+mini_side_m), ...
        (round(0.5*size(Y_im_g,2))-mini_side_n):(round(0.5*size(Y_im_g,2))+mini_side_n));
    if win
        psf_mat(((i_psf-1)*psf_side_m+1):(i_psf*psf_side_m),1:psf_side_n) = ...
            han.*psf_mat(((i_psf-1)*psf_side_m+1):(i_psf*psf_side_m),1:psf_side_n);
    end
end
% Normalization of psf_mat (for 1 PSF only):
psf_mat_norm = 1/sum(psf_mat(:))*psf_mat;
psf_struct = {psf_mat_norm,psf_size};
%
% Plotting the Kernel:
figure
imagesc(psf_mat_norm);
colormap gray;
colorbar;
title(['Kernel ',int2str(psf_side_m),'x',int2str(psf_side_n),', win = ',int2str(win)]);

%% De-convolution (No Loop):
Yrfmax = zeros(size(X1gray0));
Yrfavg = zeros(size(X1gray0));
Xmax = zeros(size(X1gray0));
Xavg = zeros(size(X1gray0));
if ~exist('numofframes','var')
    numofframes = 10; % Default.
end
Yrf_mat = zeros([size(X1gray0),numofframes]);
Yrec_mat = zeros([size(X1gray0),numofframes]);
% Parralel enhancement ("for" -> "parfor"):
frame4  = zeros([size(X1gray0),size(vid,3),numofframes]);
frame  = zeros([size(X1gray0),numofframes]);
X1gray = zeros([size(X1gray0),numofframes]);
% Stops previous parallel pool sessions:
delete(gcp('nocreate'));
% Sets the parallel pool:
% Default: parpool(4, 'IdleTimeout', 30);
%pool = parpool(4, 'IdleTimeout', Inf); 
parpool('IdleTimeout', Inf);
% Starts the parallel pool:
parfor frame_i = 1:numofframes % parfor frame_i = 1:numofframes
    % Original image normalization (from uint8[0,255] to double(0,1)):
    if class(vid) == "double"
        if size(vid,4) > 0
            frame4(:,:,:,frame_i) = vid(:,:,:,frame_i);
            if (size(frame4(:,:,:,frame_i),3) > 1)&&(max(frame4(:,:,:,frame_i),[],'all') > 1)
                X1gray(:,:,frame_i) = im2double(rgb2gray(uint8(frame4(:,:,:,frame_i)))); % Returns double [0,1].
            elseif size(frame4(:,:,:,frame_i),3) > 1
                X1gray(:,:,frame_i) = im2double(rgb2gray(frame4(:,:,:,frame_i))); % Returns double [0,1].
            elseif max(frame4(:,:,:,frame_i),[],'all') > 1
                X1gray(:,:,frame_i) = squeeze(im2double(uint8(frame4(:,:,:,frame_i)))); % Returns double [0,1].
            else
                X1gray(:,:,frame_i) = squeeze(frame4(:,:,:,frame_i)); % Returns double [0,1].
            end
        else
            frame(:,:,frame_i) = vid(:,:,frame_i);
            X1gray(:,:,frame_i) = im2double(frame(:,:,frame_i)); % Returns double [0,1].
        end
    elseif class(vid) == "VideoReader"
        frame(:,:,frame_i) = read(vid,frame_i);
        if size(frame(:,:,frame_i),3) == 1
            X1gray(:,:,frame_i) = im2double(frame(:,:,frame_i)); % Returns double [0,1].
        else
            X1gray(:,:,frame_i) = im2double(rgb2gray(frame(:,:,frame_i))); % Returns double [0,1].
        end
    else
        error('Wrong class input.\n');
    end
    
    % Instead of "im2double":
    % X1gray(:,:,frame_i) = 1/255*X1gray(:,:,frame_i);
    
    % The deconvolution algorithm:
    %[Yrf1(:,:,frame_i),Yrec1(:,:,frame_i)] = AlgFunc5(X1gray(:,:,frame_i),psf_struct,state,fil,sat,param_struct);
    %Yrf_mat(:,:,frame_i) = Yrf1(:,:,frame_i);
    %Yrec_mat(:,:,frame_i) = Yrec1(:,:,frame_i);
    [Yrf_mat(:,:,frame_i),Yrec_mat(:,:,frame_i)] = AlgFunc5_2(X1gray(:,:,frame_i),psf_struct,state,fil,sat,param_struct);
    
    % Maximal value:
    %Yrfmax = max(Yrfmax,Yrf1);
    %Xmax = max(Xmax,X1gray(:,:,frame_i));
    
    % Average value:
    %Yrfavg = Yrfavg + Yrf1;
    %Xavg = Xavg + X1gray(:,:,frame_i);
    
    % Plot Frame i:
    %{
    % RF
    figure
    image(Yrf1);
    colormap gray;
    title(['Yrf1 ',ver,', \sigma_m = ',num2str(m_sig2),', \sigma_n = ',num2str(n_sig2),', frame ',int2str(frame_i)]);
    saveas(gcf,['Yrf',int2str(frame_i),'_',ver,'.jpg']);
    % RF scaled
    figure
    imagesc(Yrf1);
    colormap gray;
    title(['Yrf1 ',ver,', \sigma_m = ',num2str(m_sig2),', \sigma_n = ',num2str(n_sig2),', scaled, frame ',int2str(frame_i)]);
    saveas(gcf,['Yrf',int2str(frame_i),'_',ver,'_sc.jpg']);
    %}
    %{
    % Recon
    figure
    image(Yrec1);
    colormap gray;
    title(['Yrec1 ',ver,', \sigma_m = ',num2str(m_sig2),', \sigma_n = ',num2str(n_sig2),', frame ',int2str(frame_i)]);
    saveas(gcf,['Yrec',int2str(frame_i),'_',ver,'.jpg']);
    % Recon scaled
    figure
    imagesc(Yrec1);
    colormap gray;
    title(['Yrec1 ',ver,', \sigma_m = ',num2str(m_sig2),', \sigma_n = ',num2str(n_sig2),', scaled, frame ',int2str(frame_i)]);
    saveas(gcf,['Yrec',int2str(frame_i),'_',ver,'_sc.jpg']);
    %}
    
    % Plot Frame i:
    %{
    % Original:
    figure
    image(X1gray(:,:,frame_i));
    colormap gray;
    title(['X1gray, frame ',int2str(frame_i)]);
    saveas(gcf,['X',int2str(frame_i),'.jpg']);
    %
    % Original scaled:
    figure
    imagesc(X1gray(:,:,frame_i));
    colormap gray;
    title(['X1gray, scaled, frame ',int2str(frame_i)]);
    saveas(gcf,['X',int2str(frame_i),'_sc.jpg']);
    %}
    
end
% Stops the parallel pool:
%pool.IdleTimeout = 30;
delete(gcp('nocreate'));

% Maximal value:
Yrfmax = max(Yrf_mat,[],3);
Xmax = max(X1gray,[],3);

% Average value:
Yrfavg = mean(Yrf_mat,3);
Xavg = mean(X1gray,3);

%Yrfavg = Yrfavg/numofframes;
%Xavg = Xavg/numofframes;

% Plot All Frames:
% Total RF sum
%{
figure
image(Yrfmax);
colormap gray;
title(['Yrfmax ',ver,', \sigma_m = ',num2str(m_sig2),', \sigma_n = ',num2str(n_sig2)]);
saveas(gcf,['Yrfmax_',ver,'.jpg']);
%}
% Total RF sum scaled
%{
figure
imagesc(Yrfmax);
colormap gray;
title(['Yrfmax ',ver,', \sigma_m = ',num2str(m_sig2),', \sigma_n = ',num2str(n_sig2),', scaled']);
saveas(gcf,['Yrfmax_',ver,'_sc.jpg']);
%}
% Total RF average
%{
figure
image(Yrfavg);
colormap gray;
title(['Yrfavg ',ver,', \sigma_m = ',num2str(m_sig2),', \sigma_n = ',num2str(n_sig2)]);
saveas(gcf,['Yrfavg_',ver,'.jpg']);
%}
% Total RF average scaled
%{
figure
imagesc(Yrfavg);
colormap gray;
title(['Yrfavg ',ver,', \sigma_m = ',num2str(m_sig2),', \sigma_n = ',num2str(n_sig2),', scaled']);
saveas(gcf,['Yrfavg_',ver,'_sc.jpg']);
%}
% Maximum Intensity Persistence
%{
figure
image(Xmax);
colormap gray;
title('Xmax - MIP');
saveas(gcf,'Xmax.jpg');
%}
% Maximum Intensity Persistence scaled
%{
figure
imagesc(Xmax);
colormap gray;
title('Xmax - MIP, scaled');
saveas(gcf,'Xmax_sc.jpg');
%}
% Total RF average
%{
figure
image(Xavg);
colormap gray;
title('Xavg - Mean');
saveas(gcf,'Xavg.jpg');
%}
% Total RF average scaled
%{
figure
imagesc(Xavg);
colormap gray;
title('Xavg - Mean, scaled');
saveas(gcf,'Xavg_sc.jpg');
%}
end