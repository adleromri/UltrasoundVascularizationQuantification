%% Image Registration:
% Registering the images according to the base frame.
%
% Syntax:
% [shortvidreg,shortvidreg_2] = PanGUI_Reg(shortvidout_pool,shortvidout_2_pool,ind_reg,N,AFS,PL)
%
% Input:
% shortvidout_pool - Original main video (B-Mode) double[0,255].
% shortvidout_2_pool - Original second video (CEUS) double[0,255].
% ind_reg - The index of the base frame.
% N - Number of iterations. (Optional)
% AFS - Accumulated Field Smoothing each iteration. (Optional)
% PL - Number of Pyramid levels. (Optional)
%
% Output:
% shortvidreg - Registration of the main video (B-Mode) double[0,255].
% shortvidreg_2 - Registration of the second video (CEUS) double[0,255].

function [shortvidreg,shortvidreg_2] = PanGUI_Reg(shortvidout_pool,shortvidout_2_pool,ind_reg,N,AFS,PL)
% Registration process:
% Parameters:
% Linear interpolation is used, as in Hanna's case.
if ~exist('N','var')
    N = 3; % default: "N = 100;", 3 in Hanna's case.
end
if ~exist('AFS','var')
    AFS = 2.00; % default: "AFS = 1.0;", 1.00 in Hanna's case?
end
if ~exist('PL','var')
    PL = 4; % default: "PL = 3;", 4 in Hanna's case.
end

% Registration 1:
% bmode:
shortvidreg_pos = shortvidout_pool(:,:,1,ind_reg:end);
L_pos = size(shortvidreg_pos,4);
shortvidreg_neg = flip(shortvidout_pool(:,:,1,1:ind_reg),4);
L_neg = size(shortvidreg_neg,4);
% ceus:
shortvidreg_pos_2 = shortvidout_2_pool(:,:,1,ind_reg:end);
shortvidreg_neg_2 = flip(shortvidout_2_pool(:,:,1,1:ind_reg),4);
%shortvidreg_pos2 = shortvid(:,:,1,ind_reg:end); % see if it works
%shortvidreg_neg2 = flip(shortvid(:,:,1,1:ind_reg),4); % see if it works

% Saves memory space in workspace:
%clear shortvidout_pool shortvidout_2_pool

for i = 1:L_pos-1
    % bmode:
    moving_i = shortvidreg_pos(:,:,1,i+1);
    %fixed_i = shortvidreg_pos(:,:,1,i);
    fixed_i = shortvidreg_pos(:,:,1,1);
    [D_i,moving_reg_i] = imregdemons(moving_i,fixed_i,N,'AccumulatedFieldSmoothing',AFS,'PyramidLevels',PL);
    shortvidreg_pos(:,:,1,i+1) = moving_reg_i;
    % ceus:
    moving_i_2 = shortvidreg_pos_2(:,:,1,i+1);
    %fixed_i_2 = shortvidreg_pos_2(:,:,1,i);
    shortvidreg_pos_2(:,:,1,i+1) = imwarp(moving_i_2,D_i);
end
for i = 1:L_neg-1
    % bmode:
    moving_i = shortvidreg_neg(:,:,1,i+1);
    %fixed_i = shortvidreg_neg(:,:,1,i);
    fixed_i = shortvidreg_neg(:,:,1,1);
    [D_i,moving_reg_i] = imregdemons(moving_i,fixed_i,N,'AccumulatedFieldSmoothing',AFS,'PyramidLevels',PL);
    shortvidreg_neg(:,:,1,i+1) = moving_reg_i;
    % ceus:
    moving_i_2 = shortvidreg_neg_2(:,:,1,i+1);
    %fixed_i_2 = shortvidreg_neg_2(:,:,1,i);
    shortvidreg_neg_2(:,:,1,i+1) = imwarp(moving_i_2,D_i);
end

% Combine the 2 parts:
% bmode:
shortvidreg = cat(4,flip(shortvidreg_neg(:,:,1,2:end),4),shortvidreg_pos);
% ceus:
shortvidreg_2 = cat(4,flip(shortvidreg_neg_2(:,:,1,2:end),4),shortvidreg_pos_2);

% Saves memory space in workspace:
clear shortvidreg_neg shortvidreg_pos
clear shortvidreg_neg_2 shortvidreg_pos_2

% Updating video to uint8[0,255] values:
% Intensity correction:
% bmode:
%shortvidreg = uint8(cat(3,shortvidreg,cat(3,shortvidreg,shortvidreg)));
%shortvidout = uint8(cat(3,shortvidreg,cat(3,shortvidreg,shortvidreg)));
shortvidreg = double(cat(3,shortvidreg,cat(3,shortvidreg,shortvidreg)));
% ceus:
%shortvidreg_2 = uint8(cat(3,shortvidreg_2,cat(3,shortvidreg_2,shortvidreg_2)));
shortvidreg_2 = double(cat(3,shortvidreg_2,cat(3,shortvidreg_2,shortvidreg_2)));
%{
if strcmp(gamma,'I')
    shortvidout_2_int = uint8(cat(3,shortvidreg_2,cat(3,shortvidreg_2,shortvidreg_2)));
    shortvidout_2 = shortvidout_2_int;
else
    shortvidout_2 = uint8(cat(3,shortvidreg_2,cat(3,shortvidreg_2,shortvidreg_2)));
end
%}

% Saves the video:
%{
% bmode:
%v = VideoWriter([respar,'Registered - ',name],'Grayscale AVI'); % Saturated video
v = VideoWriter([respar,'Registered_B - ',name],'Uncompressed AVI'); % Saturated video
v.FrameRate = fs;
open(v);
writeVideo(v,shortvidout);
close(v);
% ceus:
v = VideoWriter([respar,'Registered_CEUS - ',name],'Uncompressed AVI'); % Saturated video
v.FrameRate = fs;
open(v);
writeVideo(v,shortvidout_2);
close(v);

% Saves Original videos:
save([respar,'Registered_B - ',name,'.mat'],'shortvidreg');
save([respar,'Registered_CEUS - ',name,'.mat'],'shortvidreg_2');
%}

% Saves memory space in workspace:
%clear shortvidreg shortvidreg_2
end