%% Cross-correlation (frames pooling):
% Finding the frames which are best correlated to the base frame.
%
% Syntax:
% [shortvidout_pool,val2,ind2] = PanGUI_CC(shortvidout,frame,corr_perc,corr_thr)
%
% Input:
% shortvidout - Original main video (B-Mode) double[0,255].
% frame - The base frame for correlation compare.
% corr_perc - Correlation top percentage. (Optional)
% corr_thr - Correlation match threshold. (Optional)
%
% Output:
% shortvidout_pool - Selected frames (chronically ordered) of the main video (B-Mode) double[0,255].
% val2 - Frames pooling arrangement by frames.
% ind2 - Frames pooling arrangement by indeces.

function [shortvidout_pool,val2,ind2] = PanGUI_CC(shortvidout,frame,corr_perc,corr_thr)
% Finds cross-correlation for frame match:
if ~exist('corr_perc','var')
    corr_perc = 0.1; % Correlation top percentage.
end
if ~exist('corr_thr','var')
    corr_thr = 0.66; % Correlation match threshold.
end
all_cor = zeros(size(shortvidout,4),1);
for i = 1:size(shortvidout,4)
    all_cor_temp = normxcorr2(shortvidout(:,:,1,frame),shortvidout(:,:,1,i));
    all_cor(i) = max(all_cor_temp(:));
end
[val,ind] = sort(all_cor,'descend');
ind_top = ind(1:round(corr_perc*size(shortvidout,4)));
val_top = val(1:round(corr_perc*size(shortvidout,4)));
% Chronical (temporal) ordering:
[val2,ind2] = sort(ind_top(val_top>corr_thr));
shortvidout_pool = shortvidout(:,:,:,val2); % Chronically ordered.
%shortvidout_pool = shortvidout(:,:,:,ind_top(val_top>corr_thr)); % Ordered by matching.
end