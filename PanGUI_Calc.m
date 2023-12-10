%% Plot Calculation:
% Calculating the data for ploting.
%
% Syntax:
%{
[frames,...
    res_mean1,res_std1,res_mean1_CEUS,res_std1_CEUS,...
    res_mean2,res_std2,res_mean2_CEUS,res_std2_CEUS,...
    res_mean3,res_std3,res_mean3_CEUS,res_std3_CEUS,...
    res_mean4,res_std4,res_mean4_CEUS,res_std4_CEUS,...
    img_all,img_all_CEUS,...
    res_mean_all,res_std_all,res_mean_all_CEUS,res_std_all_CEUS,...
    meanmean1,meanstd1,meanmean1_CEUS,meanstd1_CEUS,...
    meanmean2,meanstd2,meanmean2_CEUS,meanstd2_CEUS,...
    meanmean3,meanstd3,meanmean3_CEUS,meanstd3_CEUS,...
    meanmean4,meanstd4,meanmean4_CEUS,meanstd4_CEUS,...
    mean_mean_all,mean_std_all,mean_mean_all_CEUS,mean_std_all_CEUS] = PanGUI_Calc(shortvidreg,shortvidreg_2,L,img_L)
%}
%
% Input:
% shortvidreg - Original main video (B-Mode).
% shortvidreg_2 - Original second video (CEUS).
% L - Video's length (shortvidreg).
% img_L - Labeled segmented areas.
%
% Output:
% ...

function frames = PanGUI_Calc(shortvidreg,shortvidreg_2,L,img_L,description)
%{
[frames,...
    res_mean1,res_std1,res_mean1_CEUS,res_std1_CEUS,...
    res_mean2,res_std2,res_mean2_CEUS,res_std2_CEUS,...
    res_mean3,res_std3,res_mean3_CEUS,res_std3_CEUS,...
    res_mean4,res_std4,res_mean4_CEUS,res_std4_CEUS,...
    img_all,img_all_CEUS,...
    res_mean_all,res_std_all,res_mean_all_CEUS,res_std_all_CEUS,...
    meanmean1,meanstd1,meanmean1_CEUS,meanstd1_CEUS,...
    meanmean2,meanstd2,meanmean2_CEUS,meanstd2_CEUS,...
    meanmean3,meanstd3,meanmean3_CEUS,meanstd3_CEUS,...
    meanmean4,meanstd4,meanmean4_CEUS,meanstd4_CEUS,...
    mean_mean_all,mean_std_all,mean_mean_all_CEUS,mean_std_all_CEUS] = PanGUI_Calc(shortvidreg,shortvidreg_2,L,img_L)
%}
% Define matrices:
res_mean1 = zeros(L,1);
res_std1 = zeros(L,1);
res_mean1_CEUS = zeros(L,1);
res_std1_CEUS = zeros(L,1);
res_mean2 = zeros(L,1);
res_std2 = zeros(L,1);
res_mean2_CEUS = zeros(L,1);
res_std2_CEUS = zeros(L,1);
res_mean3 = zeros(L,1);
res_std3 = zeros(L,1);
res_mean3_CEUS = zeros(L,1);
res_std3_CEUS = zeros(L,1);
res_mean4 = zeros(L,1);
res_std4 = zeros(L,1);
res_mean4_CEUS = zeros(L,1);
res_std4_CEUS = zeros(L,1);
% [104:495,1:540],[1:103,79:468]
res_mean_all = zeros(L,1);
res_std_all = zeros(L,1);
res_mean_all_CEUS = zeros(L,1);
res_std_all_CEUS = zeros(L,1);

% Results calculation:
for i = 1:L
    % Reads the frame:
    img = shortvidreg(:,:,1,i);
    img_CEUS = shortvidreg_2(:,:,1,i);
    % Converting values from uint8[0,255] to double[0,1].
    img_res = 1/255*double(img(:,:,1));
    img_res_CEUS = 1/255*double(img_CEUS(:,:,1));
    
    % Analyzes:
    res_mean1(i) = mean(img_res(img_L == 1));
    res_std1(i) = std(img_res(img_L == 1));
    res_mean1_CEUS(i) = mean(img_res_CEUS(img_L == 1));
    res_std1_CEUS(i) = std(img_res_CEUS(img_L == 1));
    
    res_mean2(i) = mean(img_res(img_L == 2));
    res_std2(i) = std(img_res(img_L == 2));
    res_mean2_CEUS(i) = mean(img_res_CEUS(img_L == 2));
    res_std2_CEUS(i) = std(img_res_CEUS(img_L == 2));
    
    res_mean3(i) = mean(img_res(img_L == 3));
    res_std3(i) = std(img_res(img_L == 3));
    res_mean3_CEUS(i) = mean(img_res_CEUS(img_L == 3));
    res_std3_CEUS(i) = std(img_res_CEUS(img_L == 3));
    
    res_mean4(i) = mean(img_res(img_L == 4));
    res_std4(i) = std(img_res(img_L == 4));
    res_mean4_CEUS(i) = mean(img_res_CEUS(img_L == 4));
    res_std4_CEUS(i) = std(img_res_CEUS(img_L == 4));
    
    % [104:495,1:540],[1:103,79:468]
    %img_all = [reshape(img_res(104:495,1:540),(495-104+1)*(540-1+1),1);reshape(img_res(1:103,79:468),(103-1+1)*(468-79+1),1)];
    %img_all_CEUS = [reshape(img_res_CEUS(104:495,1:540),(495-104+1)*(540-1+1),1);reshape(img_res_CEUS(1:103,79:468),(103-1+1)*(468-79+1),1)];
    % No upper part (near field):
    img_all = reshape(img_res(104:495,1:540),(495-104+1)*(540-1+1),1);
    img_all_CEUS = reshape(img_res_CEUS(104:495,1:540),(495-104+1)*(540-1+1),1);
    res_mean_all(i) = mean(img_all);
    res_std_all(i) = std(img_all);
    res_mean_all_CEUS(i) = mean(img_all_CEUS);
    res_std_all_CEUS(i) = std(img_all_CEUS);
end

% Prepares data for plotting:
frames = 1:L;
meanmean1 = mean(res_mean1(1:L));
meanstd1 = std(res_mean1(1:L));
meanmean1_CEUS = mean(res_mean1_CEUS(1:L));
meanstd1_CEUS = std(res_mean1_CEUS(1:L));
meanmean2 = mean(res_mean2(1:L));
meanstd2 = std(res_mean2(1:L));
meanmean2_CEUS = mean(res_mean2_CEUS(1:L));
meanstd2_CEUS = std(res_mean2_CEUS(1:L));
meanmean3 = mean(res_mean3(1:L));
meanstd3 = std(res_mean3(1:L));
meanmean3_CEUS = mean(res_mean3_CEUS(1:L));
meanstd3_CEUS = std(res_mean3_CEUS(1:L));
meanmean4 = mean(res_mean4(1:L));
meanstd4 = std(res_mean4(1:L));
meanmean4_CEUS = mean(res_mean4_CEUS(1:L));
meanstd4_CEUS = std(res_mean4_CEUS(1:L));
% [104:495,1:540],[1:103,79:468]
mean_mean_all = mean(res_mean_all(1:L));
mean_std_all = std(res_mean_all(1:L));
mean_mean_all_CEUS = mean(res_mean_all_CEUS(1:L));
mean_std_all_CEUS = std(res_mean_all_CEUS(1:L));

% Data smoothing:
% LPF:
%{
x = res_mean1;
wpass = 0.5;
y = lowpass(x,wpass);
y2 = y+mean(x)-mean(y);
%}
% Polynom:
%{
n = 6;
x = (frames).';
p1 = polyfit(x,res_mean1,n);
p2 = polyfit(x,res_mean2,n);
p3 = polyfit(x,res_mean3,n);
p4 = polyfit(x,res_mean4,n);
p5 = polyfit(x,res_mean_all,n);
pceus1 = polyfit(x,res_mean1_CEUS,n);
pceus2 = polyfit(x,res_mean2_CEUS,n);
pceus3 = polyfit(x,res_mean3_CEUS,n);
pceus4 = polyfit(x,res_mean4_CEUS,n);
pceus5 = polyfit(x,res_mean_all_CEUS,n);
res_mean1_pol = polyval(p1,x);
res_mean2_pol = polyval(p2,x);
res_mean3_pol = polyval(p3,x);
res_mean4_pol = polyval(p4,x);
res_mean_all_pol = polyval(p5,x);
res_mean1_CEUS_pol = polyval(pceus1,x);
res_mean2_CEUS_pol = polyval(pceus2,x);
res_mean3_CEUS_pol = polyval(pceus3,x);
res_mean4_CEUS_pol = polyval(pceus4,x);
res_mean_all_CEUS_pol = polyval(pceus5,x);
%}
% Average:
n = 5;
ker = 1/n*ones(n,1);
res_mean1_pol = conv(res_mean1,ker,'same');
res_mean2_pol = conv(res_mean2,ker,'same');
res_mean3_pol = conv(res_mean3,ker,'same');
res_mean4_pol = conv(res_mean4,ker,'same');
res_mean_all_pol = conv(res_mean_all,ker,'same');
res_mean1_CEUS_pol = conv(res_mean1_CEUS,ker,'same');
res_mean2_CEUS_pol = conv(res_mean2_CEUS,ker,'same');
res_mean3_CEUS_pol = conv(res_mean3_CEUS,ker,'same');
res_mean4_CEUS_pol = conv(res_mean4_CEUS,ker,'same');
res_mean_all_CEUS_pol = conv(res_mean_all_CEUS,ker,'same');
% Fixing edges:
res_mean1_pol = [res_mean1(1:n);res_mean1_pol(n+1:end-n);res_mean1(end-n+1:end)];
res_mean2_pol = [res_mean2(1:n);res_mean2_pol(n+1:end-n);res_mean2(end-n+1:end)];
res_mean3_pol = [res_mean3(1:n);res_mean3_pol(n+1:end-n);res_mean3(end-n+1:end)];
res_mean4_pol = [res_mean4(1:n);res_mean4_pol(n+1:end-n);res_mean4(end-n+1:end)];
res_mean_all_pol = [res_mean_all(1:n);res_mean_all_pol(n+1:end-n);res_mean_all(end-n+1:end)];
res_mean1_CEUS_pol = [res_mean1_CEUS(1:n);res_mean1_CEUS_pol(n+1:end-n);res_mean1_CEUS(end-n+1:end)];
res_mean2_CEUS_pol = [res_mean2_CEUS(1:n);res_mean2_CEUS_pol(n+1:end-n);res_mean2_CEUS(end-n+1:end)];
res_mean3_CEUS_pol = [res_mean3_CEUS(1:n);res_mean3_CEUS_pol(n+1:end-n);res_mean3_CEUS(end-n+1:end)];
res_mean4_CEUS_pol = [res_mean4_CEUS(1:n);res_mean4_CEUS_pol(n+1:end-n);res_mean4_CEUS(end-n+1:end)];
res_mean_all_CEUS_pol = [res_mean_all_CEUS(1:n);res_mean_all_CEUS_pol(n+1:end-n);res_mean_all_CEUS(end-n+1:end)];

% Histograms:
% Reads the videos:
% Converting values from uint8[0,255] to double[0,255].
img_res = double(shortvidreg(:,:,1,:));
img_res_CEUS = double(shortvidreg_2(:,:,1,:));
img_all = reshape(img_res(104:495,1:540,1,:),(495-104+1)*(540-1+1)*L,1,1,1);
img_all_CEUS = reshape(img_res_CEUS(104:495,1:540,1,:),(495-104+1)*(540-1+1)*L,1,1,1);
res_hist1 = img_res(img_L == 1);
res_hist2 = img_res(img_L == 2);
res_hist3 = img_res(img_L == 3);
res_hist4 = img_res(img_L == 4);
res_hist_all = img_all;
res_hist1_CEUS = img_res_CEUS(img_L == 1);
res_hist2_CEUS = img_res_CEUS(img_L == 2);
res_hist3_CEUS = img_res_CEUS(img_L == 3);
res_hist4_CEUS = img_res_CEUS(img_L == 4);
res_hist_all_CEUS = img_all_CEUS;



% Saves the data:
save(['.\Results - temp\',description,'_data.mat'],'img_all','img_all_CEUS',...
    'res_mean_all','res_std_all','res_mean_all_CEUS','res_std_all_CEUS',...
    'mean_mean_all','mean_std_all','mean_mean_all_CEUS','mean_std_all_CEUS',...
    'res_mean1','res_std1','res_mean1_CEUS','res_std1_CEUS',...
    'res_mean2','res_std2','res_mean2_CEUS','res_std2_CEUS',...
    'res_mean3','res_std3','res_mean3_CEUS','res_std3_CEUS',...
    'res_mean4','res_std4','res_mean4_CEUS','res_std4_CEUS',...
    'meanmean1','meanstd1','meanmean1_CEUS','meanstd1_CEUS',...
    'meanmean2','meanstd2','meanmean2_CEUS','meanstd2_CEUS',...
    'meanmean3','meanstd3','meanmean3_CEUS','meanstd3_CEUS',...
    'meanmean4','meanstd4','meanmean4_CEUS','meanstd4_CEUS',...
    'res_mean1_pol','res_mean2_pol','res_mean3_pol','res_mean4_pol','res_mean_all_pol',...
    'res_mean1_CEUS_pol','res_mean2_CEUS_pol','res_mean3_CEUS_pol','res_mean4_CEUS_pol','res_mean_all_CEUS_pol',...
    'res_hist1','res_hist2','res_hist3','res_hist4','res_hist_all',...
    'res_hist1_CEUS','res_hist2_CEUS','res_hist3_CEUS','res_hist4_CEUS','res_hist_all_CEUS');
end