%% Kernel Check
% Finds the biggest Kernel for deconvolution that still fits the given image.
%
% Syntax:
% kerfit = ker_check(img,param_struct)
%
% Input:
% img - the tested image (grayscale).
% param_struct - parameters structure, including: max_dif, max_loop,
% change_delay, type, nn.
%
% Output:
% kerfit - the best fit deconvolution kernel.
% sigma_m - the axial (y-axis) sigma of the gaussian kernel.
% sigma_n - the lateral (x-axis) sigma of the gaussian kernel.

function [kerfit,sigma_m,sigma_n] = ker_check_2(img,param_struct)
%% Set the parameters:
% Kernel minimal sigma value:
val_min = 0.5; % 0.5
% Kernel maximal sigma value:
val_max = 20.0; % 1.5
% Kernel sigma value resolution:
val_res = 0.1; % 0.1
% Maximal Absolute Difference limit:
MAD_thr = 0.1; % 0.5
% Average Absolute Difference limit:
AAD_thr = 0.01; % 0.05
% Minimal peak value to treat:
peak_min = 0.1; % 0.1

%% Calculations:
val = val_min:val_res:val_max;
L = length(val);
val_avg = zeros(L,1);
MAD = zeros(L,1);
AAD = zeros(L,1);
% Parralel enhancement ("for" -> "parfor"):
parfor i = 0:(L-1)
    sigma_m = val_min+i*val_res;
    sigma_n = val_min+i*val_res;
    [val_avg(i+1),MAD(i+1),AAD(i+1)] = val_check(img,sigma_m,sigma_n,peak_min,param_struct);
end
% Removes wide Kernels which reconstruct poorly:
for i = 1:L
    if (MAD(i) < MAD_thr)&&(AAD(i) < AAD_thr)
        val_avg(i) = val_avg(i);
    else
        val_avg(i) = max(val_avg);
    end
end

%% The chosen kernel:
% [val,ind] = min(val_avg);
[~,val_fit_ind] = min(val_avg);
i = val_fit_ind - 1;
sigma_m = val_min+i*val_res;
sigma_n = val_min+i*val_res;

% Gaussian Kernel (5x5):
M = 7; % Default: 7 or 5.
N = 7; % Default: 7 or 5.
ker = zeros(M,N);
max_gauss = 1;
mu_m = 0.5*(M+1);
mu_n = 0.5*(N+1);
for m = 1:M
    for n = 1:N
        ker(m,n) = max_gauss*exp(-0.5*([m;n]-[mu_m;mu_n]).'*([sigma_m 0;0 sigma_n]\([m;n]-[mu_m;mu_n])));
    end
end
% psf_mat = ker;
% psf_size = [M,N];
% psf_struct = {psf_mat,psf_size};

kerfit = ker;

fprintf(['The chosen Kernel has: sigma_m = ',num2str(sigma_m),', sigma_n = ',num2str(sigma_n),'\n']);

end


%% Checks the value of the surrounding reflectivity function peaks:
% Best "val_avg" would be: "0", and worst: "1".
% Assumes double pixel values: [0,1] or [0,255].
function [val_avg,MAD,AAD] = val_check(img,sigma_m,sigma_n,peak_min,param_struct)
% Gaussian Kernel (3x3):
M = 7; % Default: 7 or 5.
N = 7; % Default: 7 or 5.
ker = zeros(M,N);
max_gauss = 1;
mu_m = 0.5*(M+1);
mu_n = 0.5*(N+1);
for m = 1:M
    for n = 1:N
        ker(m,n) = max_gauss*exp(-0.5*([m;n]-[mu_m;mu_n]).'*([sigma_m 0;0 sigma_n]\([m;n]-[mu_m;mu_n])));
    end
end
psf_mat = ker;
psf_size = [M,N];
psf_struct = {psf_mat,psf_size};

% Need to set "state" (initial state: 'thr' - soft threshold, 'zero').
state = 'zero'; % Default 'thr'.
% Need to set "ver" (version: 'HalfCons','DecCons','SevCons','MaxCons','MoreCons','RegCons','DifVal').
% ver = 'RegCons'; % Default. Seems to bring the best result.
% Need to set "fil" (filtering: 'on','off').
fil = 'on'; % Default 'on'.
% Need to set "sat" (saturaion: 'on','off').
sat = 'off'; % Default.

[Yrf,Yrec] = AlgFunc5(img,psf_struct,state,fil,sat,param_struct);

% 4-neighbors filter:
% peakfilter = [0 -1 0;-1 4 -1;0 -1 0];
% 8-neighbors filter:
peakfilter = [-1 -1 -1;-1 8 -1;-1 -1 -1];

sources = conv2(Yrf,peakfilter,'same');
sources_bin = double((sources > 0)&(Yrf > peak_min));

num = sum(sources_bin(:));
val = zeros(num,1);
i = 1;
for m = 1:size(sources_bin,1)
    for n = 1:size(sources_bin,2)
        if sources_bin(m,n) == 1
            if (m > 1)&&(n > 1)&&(m < size(sources_bin,1))&&(n < size(sources_bin,2))
                val(i) = mean([Yrf(m+1,n),Yrf(m-1,n),Yrf(m,n+1),Yrf(m,n-1),...
                    Yrf(m+1,n+1),Yrf(m-1,n+1),Yrf(m+1,n-1),Yrf(m-1,n-1)])/Yrf(m,n);
            else
                val(i) = [];
                i = i - 1;
            end
            i = i + 1;
        end
    end
end

val_avg = mean(val);

MAD = max(abs(img(:)-Yrec(:))); % Maximum
% RMS = rms(abs(img(:)-Yrec(:))); % RMS
% MAD2 = median(abs(img(:)-Yrec(:))); % Median
AAD = mean(abs(img(:)-Yrec(:))); % Average

end