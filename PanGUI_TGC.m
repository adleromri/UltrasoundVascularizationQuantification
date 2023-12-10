%% Gamma correction (TGC - Time Gain Compensation):
% Increasing the image intensity with Time Gain Compensation based on Gamma correction.
%
% Syntax:
% shortvidout_2_int = PanGUI_TGC(shortvid,shortvid_2,frame)
%
% Input:
% shortvid - Original first video (B-Mode) double[0,255].
% shortvid_2 - Original second video (CEUS) double[0,255].
% frame - The selected frame.
%
% Output:
% shortvidout_2_int - Intensified second video (CEUS) uint8[0,255].

function shortvidout_2_int = PanGUI_TGC(shortvid,shortvid_2,frame)
% Origin (transducer location) and edge of the transmission:
m_edge = size(shortvid_2,1);
n_edge = size(shortvid_2,2);
L = size(shortvid_2,4);
m_origin = -n_edge/4; % m_origin = 1;
n_origin = n_edge/2 + 0.5;
eps = round(1/10*sqrt((m_edge-m_origin)^2 + (n_edge-n_origin)^2));
dist_max = sqrt((m_edge-m_origin)^2 + (n_edge-n_origin)^2) + eps;

% Convert from int[0,255] to double[0,1]:
shortvid_double = 1/255*double(shortvid);
shortvid_2_double = 1/255*double(shortvid_2);
shortvid_2_gamma = zeros(size(shortvid_2_double));

% Histogram evaluation:
secnum = 5;
secrad = dist_max/secnum;
secvec = [0,secrad,2*secrad,3*secrad,4*secrad,5*secrad];
hist1 = [];
hist2 = [];
hist3 = [];
hist4 = [];
hist5 = [];
for m = 1:size(shortvid_2,1)
    for n = 1:size(shortvid_2,2)
        dist = sqrt((m-m_origin)^2 + (n-n_origin)^2);
        if dist>=secvec(1) && dist<=secvec(2)
            %hist1 = [hist1,shortvid_2_double(m,n,1,frame)];
            if frame > 1 && frame < L
                hist1 = [hist1,shortvid_2_double(m,n,1,frame-1:frame+1)];
            elseif frame == 1
                hist1 = [hist1,shortvid_2_double(m,n,1,frame:frame+2)];
            elseif frame == L
                hist1 = [hist1,shortvid_2_double(m,n,1,frame-2:frame)];
            end
        elseif dist>secvec(2) && dist<=secvec(3)
            %hist2 = [hist2,shortvid_2_double(m,n,1,frame)];
            if frame > 1 && frame < L
                hist2 = [hist2,shortvid_2_double(m,n,1,frame-1:frame+1)];
            elseif frame == 1
                hist2 = [hist2,shortvid_2_double(m,n,1,frame:frame+2)];
            elseif frame == L
                hist2 = [hist2,shortvid_2_double(m,n,1,frame-2:frame)];
            end
        elseif dist>secvec(3) && dist<=secvec(4)
            %hist3 = [hist3,shortvid_2_double(m,n,1,frame)];
            if frame > 1 && frame < L
                hist3 = [hist3,shortvid_2_double(m,n,1,frame-1:frame+1)];
            elseif frame == 1
                hist3 = [hist3,shortvid_2_double(m,n,1,frame:frame+2)];
            elseif frame == L
                hist3 = [hist3,shortvid_2_double(m,n,1,frame-2:frame)];
            end
        elseif dist>secvec(4) && dist<=secvec(5)
            %hist4 = [hist4,shortvid_2_double(m,n,1,frame)];
            if frame > 1 && frame < L
                hist4 = [hist4,shortvid_2_double(m,n,1,frame-1:frame+1)];
            elseif frame == 1
                hist4 = [hist4,shortvid_2_double(m,n,1,frame:frame+2)];
            elseif frame == L
                hist4 = [hist4,shortvid_2_double(m,n,1,frame-2:frame)];
            end
        elseif dist>secvec(5) && dist<=secvec(6)
            %hist5 = [hist5,shortvid_2_double(m,n,1,frame)];
            if frame > 1 && frame < L
                hist5 = [hist5,shortvid_2_double(m,n,1,frame-1:frame+1)];
            elseif frame == 1
                hist5 = [hist5,shortvid_2_double(m,n,1,frame:frame+2)];
            elseif frame == L
                hist5 = [hist5,shortvid_2_double(m,n,1,frame-2:frame)];
            end
        end
    end
end
% Median version:
%{
hist1 = hist1(hist1 > 0);
hist2 = hist2(hist2 > 0);
hist3 = hist3(hist3 > 0);
hist4 = hist4(hist4 > 0);
hist5 = hist5(hist5 > 0);
%histmax = [max(hist1),max(hist2),max(hist3),max(hist4),max(hist5)];
histmed = [median(hist1),median(hist2),median(hist3),median(hist4),median(hist5)];
%histmean = [mean(hist1),mean(hist2),mean(hist3),mean(hist4),mean(hist5)];
histgamma = log(max(histmed))./log(histmed);
%}
% Mean version:
%{
eps = 0.001;
hist1 = hist1 + eps;
hist2 = hist2 + eps;
hist3 = hist3 + eps;
hist4 = hist4 + eps;
hist5 = hist5 + eps;
%}
hist1 = hist1(hist1 > 0);
hist2 = hist2(hist2 > 0);
hist3 = hist3(hist3 > 0);
hist4 = hist4(hist4 > 0);
hist5 = hist5(hist5 > 0);
histmean = [mean(hist1),mean(hist2),mean(hist3),mean(hist4),mean(hist5)];
histmean_bmode = mean(shortvid_double(:,:,1,frame),'all');
histgamma = log(histmean_bmode)./log(histmean);
%
p = polyfit(0:1/(secnum-1):1,histgamma,secnum-1); % p(1)*x^4+p(2)*x^3+p(3)*x^2+p(4)*x+p(5)
%y = polyval(p,0:1/100:1);
%{
figure
plot(0:1/(secnum-1):1,histmax);
title('max vs. depth');
figure
plot(0:1/(secnum-1):1,histmed);
title('median vs. depth');
figure
plot(0:1/(secnum-1):1,histmean);
title('mean vs. depth');
figure
plot(0:1/(secnum-1):1,histgamma,0:1/100:1,y);
legend('gamma by histogram','gamma by polyfit');
title('gamma vs. depth');
figure
hist1_0 = histogram(hist1,0:(1/255):1,'Normalization','probability');
figure
hist2_0 = histogram(hist2,0:(1/255):1,'Normalization','probability');
figure
hist3_0 = histogram(hist3,0:(1/255):1,'Normalization','probability');
figure
hist4_0 = histogram(hist4,0:(1/255):1,'Normalization','probability');
figure
hist5_0 = histogram(hist5,0:(1/255):1,'Normalization','probability');
%}
%Gamma = abs(dist_max - dist)/dist_max;
%shortvid_2_gamma(m,n,1,:) = shortvid_2_double(m,n,1,:).^Gamma;

% Gamma:
for m = 1:size(shortvid_2,1)
    for n = 1:size(shortvid_2,2)
        dist = sqrt((m-m_origin)^2 + (n-n_origin)^2);
        % Linear TGC / Gamma correction:
        %Gamma = abs(dist_max - dist)/dist_max;
        % Adaptive TGC / Gamma correction:
        x = dist/dist_max;
        Gamma = polyval(p,x);
        shortvid_2_gamma(m,n,1,:) = shortvid_2_double(m,n,1,:).^Gamma;
    end
end

% Convert from 1 [double] to 256 [int]:
shortvid_2_gamma = 255*shortvid_2_gamma;

% Preparing for "VideoWriter":
% Intensity correction:
%shortvidout_2_int = uint8(cat(3,shortvid_2_gamma,cat(3,shortvid_2_gamma,shortvid_2_gamma))); % CEUS
shortvidout_2_int = shortvid_2_gamma; % CEUS

% Saves the video:
%{
if ~exist('.\Results - temp','dir')
    mkdir('.\Results - temp');
end
% CEUS
%v = VideoWriter([respar,'Focused_2 - ',name],'Grayscale AVI');
v = VideoWriter(['.\Results - temp\Focused_CEUS + Gamma - ',name],'Uncompressed AVI');
v.FrameRate = fs;
open(v);
writeVideo(v,shortvidout_2_int);
close(v);

% Saves Original videos:
save(['.\Results - temp\TGC_CEUS - ',name,'.mat'],'shortvid_2_gamma');
%}
end