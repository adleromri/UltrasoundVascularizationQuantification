%% Image Plot:
% Ploting the analysis results.
%
% Syntax:
%{
PanGUI_Plot(img_L,pixel_cm,L,frames,img_all,img_all_CEUS,...
    meanmean1,meanstd1,res_mean1,res_std1,...
    meanmean2,meanstd2,res_mean2,res_std2,...
    meanmean3,meanstd3,res_mean3,res_std3,...
    meanmean4,meanstd4,res_mean4,res_std4,...
    meanmean1_CEUS,meanstd1_CEUS,res_mean1_CEUS,res_std1_CEUS,...
    meanmean2_CEUS,meanstd2_CEUS,res_mean2_CEUS,res_std2_CEUS,...
    meanmean3_CEUS,meanstd3_CEUS,res_mean3_CEUS,res_std3_CEUS,...
    meanmean4_CEUS,meanstd4_CEUS,res_mean4_CEUS,res_std4_CEUS,...
    mean_mean_all,mean_std_all,res_mean_all,res_std_all,...
    mean_mean_all_CEUS,mean_std_all_CEUS,res_mean_all_CEUS,res_std_all_CEUS,...
    shortvidreg,shortvidreg_2,ind_reg,img_3,img_3_CEUS,frame)
%}
%
% Input:
% ...
%
% Output:
% Saves plots in the directory: ".\Results - temp".

function PanGUI_Plot_2(img_L,pixel_cm,L,frames,...
    shortvidreg,shortvidreg_2,ind_reg,img_3,img_3_CEUS,frame,description,description_2)
%{
PanGUI_Plot(img_L,pixel_cm,L,frames,img_all,img_all_CEUS,...
    meanmean1,meanstd1,res_mean1,res_std1,...
    meanmean2,meanstd2,res_mean2,res_std2,...
    meanmean3,meanstd3,res_mean3,res_std3,...
    meanmean4,meanstd4,res_mean4,res_std4,...
    meanmean1_CEUS,meanstd1_CEUS,res_mean1_CEUS,res_std1_CEUS,...
    meanmean2_CEUS,meanstd2_CEUS,res_mean2_CEUS,res_std2_CEUS,...
    meanmean3_CEUS,meanstd3_CEUS,res_mean3_CEUS,res_std3_CEUS,...
    meanmean4_CEUS,meanstd4_CEUS,res_mean4_CEUS,res_std4_CEUS,...
    mean_mean_all,mean_std_all,res_mean_all,res_std_all,...
    mean_mean_all_CEUS,mean_std_all_CEUS,res_mean_all_CEUS,res_std_all_CEUS,...
    shortvidreg,shortvidreg_2,ind_reg,img_3,img_3_CEUS,frame)
%}
% Loads the data:
load(['.\Results - temp\',description,'_data.mat']);

% Saves labeled data and results:
% B-Mode:
f = figure('visible','off');
hold on;
area([meanmean1*ones(L,1),meanstd1*ones(L,1)],meanmean1-meanstd1)
plot(frames,res_mean1(1:L),'LineWidth',3);
errorbar(frames,res_mean1(1:L),res_std1(1:L),'LineWidth',1.5);
title(['TIC - Mean + STD, Label 1, ',int2str(sum(img_L(img_L == 1))),...
    ' pixels (~',num2str(pixel_cm^2*sum(img_L(img_L == 1))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
ylim([0,1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_1.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
area([meanmean2*ones(L,1),meanstd2*ones(L,1)],meanmean2-meanstd2)
plot(frames,res_mean2(1:L),'LineWidth',3);
errorbar(frames,res_mean2(1:L),res_std2(1:L),'LineWidth',1.5);
title(['TIC - Mean + STD, Label 2, ',int2str(sum(1/2*img_L(img_L == 2))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/2*img_L(img_L == 2))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
ylim([0,1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_2.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
area([meanmean3*ones(L,1),meanstd3*ones(L,1)],meanmean3-meanstd3)
plot(frames,res_mean3(1:L),'LineWidth',3);
errorbar(frames,res_mean3(1:L),res_std3(1:L),'LineWidth',1.5);
title(['TIC - Mean + STD, Label 3, ',int2str(sum(1/3*img_L(img_L == 3))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/3*img_L(img_L == 3))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
ylim([0,1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_3.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
area([meanmean4*ones(L,1),meanstd4*ones(L,1)],meanmean4-meanstd4)
plot(frames,res_mean4(1:L),'LineWidth',3);
errorbar(frames,res_mean4(1:L),res_std4(1:L),'LineWidth',1.5);
title(['TIC - Mean + STD, Label 4, ',int2str(sum(1/4*img_L(img_L == 4))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/4*img_L(img_L == 4))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
ylim([0,1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_4.jpg']);
close(f);
%
% CEUS:
f = figure('visible','off');
hold on;
area([meanmean1_CEUS*ones(L,1),meanstd1_CEUS*ones(L,1)],meanmean1_CEUS-meanstd1_CEUS)
plot(frames,res_mean1_CEUS(1:L),'LineWidth',3);
errorbar(frames,res_mean1_CEUS(1:L),res_std1_CEUS(1:L),'LineWidth',1.5);
title(['TIC (CEUS) - Mean + STD, Label 1, ',int2str(sum(img_L(img_L == 1))),...
    ' pixels (~',num2str(pixel_cm^2*sum(img_L(img_L == 1))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
ylim([0,1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_1_CEUS.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
area([meanmean2_CEUS*ones(L,1),meanstd2_CEUS*ones(L,1)],meanmean2_CEUS-meanstd2_CEUS)
plot(frames,res_mean2_CEUS(1:L),'LineWidth',3);
errorbar(frames,res_mean2_CEUS(1:L),res_std2_CEUS(1:L),'LineWidth',1.5);
title(['TIC (CEUS) - Mean + STD, Label 2, ',int2str(sum(1/2*img_L(img_L == 2))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/2*img_L(img_L == 2))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
ylim([0,1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_2_CEUS.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
area([meanmean3_CEUS*ones(L,1),meanstd3_CEUS*ones(L,1)],meanmean3_CEUS-meanstd3_CEUS)
plot(frames,res_mean3_CEUS(1:L),'LineWidth',3);
errorbar(frames,res_mean3_CEUS(1:L),res_std3_CEUS(1:L),'LineWidth',1.5);
title(['TIC (CEUS) - Mean + STD, Label 3, ',int2str(sum(1/3*img_L(img_L == 3))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/3*img_L(img_L == 3))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
ylim([0,1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_3_CEUS.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
area([meanmean4_CEUS*ones(L,1),meanstd4_CEUS*ones(L,1)],meanmean4_CEUS-meanstd4_CEUS)
plot(frames,res_mean4_CEUS(1:L),'LineWidth',3);
errorbar(frames,res_mean4_CEUS(1:L),res_std4_CEUS(1:L),'LineWidth',1.5);
title(['TIC (CEUS) - Mean + STD, Label 4, ',int2str(sum(1/4*img_L(img_L == 4))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/4*img_L(img_L == 4))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
ylim([0,1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_4_CEUS.jpg']);
close(f);


% Saves labeled data and results (zoom):
% B-Mode (zoom):
f = figure('visible','off');
hold on;
area([meanmean1*ones(L,1),meanstd1*ones(L,1)],meanmean1-meanstd1)
plot(frames,res_mean1(1:L),'LineWidth',3);
errorbar(frames,res_mean1(1:L),res_std1(1:L),'LineWidth',1.5);
title(['TIC - Mean + STD, Label 1, ',int2str(sum(img_L(img_L == 1))),...
    ' pixels (~',num2str(pixel_cm^2*sum(img_L(img_L == 1))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_1_zoom.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
area([meanmean2*ones(L,1),meanstd2*ones(L,1)],meanmean2-meanstd2)
plot(frames,res_mean2(1:L),'LineWidth',3);
errorbar(frames,res_mean2(1:L),res_std2(1:L),'LineWidth',1.5);
title(['TIC - Mean + STD, Label 2, ',int2str(sum(1/2*img_L(img_L == 2))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/2*img_L(img_L == 2))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_2_zoom.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
area([meanmean3*ones(L,1),meanstd3*ones(L,1)],meanmean3-meanstd3)
plot(frames,res_mean3(1:L),'LineWidth',3);
errorbar(frames,res_mean3(1:L),res_std3(1:L),'LineWidth',1.5);
title(['TIC - Mean + STD, Label 3, ',int2str(sum(1/3*img_L(img_L == 3))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/3*img_L(img_L == 3))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_3_zoom.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
area([meanmean4*ones(L,1),meanstd4*ones(L,1)],meanmean4-meanstd4)
plot(frames,res_mean4(1:L),'LineWidth',3);
errorbar(frames,res_mean4(1:L),res_std4(1:L),'LineWidth',1.5);
title(['TIC - Mean + STD, Label 4, ',int2str(sum(1/4*img_L(img_L == 4))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/4*img_L(img_L == 4))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_4_zoom.jpg']);
close(f);
%
% CEUS (zoom):
f = figure('visible','off');
hold on;
area([meanmean1_CEUS*ones(L,1),meanstd1_CEUS*ones(L,1)],meanmean1_CEUS-meanstd1_CEUS)
plot(frames,res_mean1_CEUS(1:L),'LineWidth',3);
errorbar(frames,res_mean1_CEUS(1:L),res_std1_CEUS(1:L),'LineWidth',1.5);
title(['TIC (CEUS) - Mean + STD, Label 1, ',int2str(sum(img_L(img_L == 1))),...
    ' pixels (~',num2str(pixel_cm^2*sum(img_L(img_L == 1))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_1_CEUS_zoom.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
area([meanmean2_CEUS*ones(L,1),meanstd2_CEUS*ones(L,1)],meanmean2_CEUS-meanstd2_CEUS)
plot(frames,res_mean2_CEUS(1:L),'LineWidth',3);
errorbar(frames,res_mean2_CEUS(1:L),res_std2_CEUS(1:L),'LineWidth',1.5);
title(['TIC (CEUS) - Mean + STD, Label 2, ',int2str(sum(1/2*img_L(img_L == 2))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/2*img_L(img_L == 2))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_2_CEUS_zoom.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
area([meanmean3_CEUS*ones(L,1),meanstd3_CEUS*ones(L,1)],meanmean3_CEUS-meanstd3_CEUS)
plot(frames,res_mean3_CEUS(1:L),'LineWidth',3);
errorbar(frames,res_mean3_CEUS(1:L),res_std3_CEUS(1:L),'LineWidth',1.5);
title(['TIC (CEUS) - Mean + STD, Label 3, ',int2str(sum(1/3*img_L(img_L == 3))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/3*img_L(img_L == 3))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_3_CEUS_zoom.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
area([meanmean4_CEUS*ones(L,1),meanstd4_CEUS*ones(L,1)],meanmean4_CEUS-meanstd4_CEUS)
plot(frames,res_mean4_CEUS(1:L),'LineWidth',3);
errorbar(frames,res_mean4_CEUS(1:L),res_std4_CEUS(1:L),'LineWidth',1.5);
title(['TIC (CEUS) - Mean + STD, Label 4, ',int2str(sum(1/4*img_L(img_L == 4))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/4*img_L(img_L == 4))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_4_CEUS_zoom.jpg']);
close(f);
%
% All (zoom):
f = figure('visible','off');
hold on;
area([mean_mean_all*ones(L,1),mean_std_all*ones(L,1)],mean_mean_all-mean_std_all)
plot(frames,res_mean_all(1:L),'LineWidth',3);
errorbar(frames,res_mean_all(1:L),res_std_all(1:L),'LineWidth',1.5);
title(['TIC - Mean + STD, Label All, ',int2str(sum(img_all == img_all)),...
    ' pixels (~',num2str(pixel_cm^2*sum(img_all == img_all)),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_All_zoom.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
area([mean_mean_all_CEUS*ones(L,1),mean_std_all_CEUS*ones(L,1)],mean_mean_all_CEUS-mean_std_all_CEUS)
plot(frames,res_mean_all_CEUS(1:L),'LineWidth',3);
errorbar(frames,res_mean_all_CEUS(1:L),res_std_all_CEUS(1:L),'LineWidth',1.5);
title(['TIC (CEUS) - Mean + STD, Label All, ',int2str(sum(img_all_CEUS == img_all_CEUS)),...
    ' pixels (~',num2str(pixel_cm^2*sum(img_all_CEUS == img_all_CEUS)),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_All_CEUS_zoom.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
plot(frames,res_mean1(1:L),frames,res_mean2(1:L),frames,res_mean3(1:L),...
    frames,res_mean4(1:L),frames,res_mean_all(1:L));
title('TIC - Mean All Labels');
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
%ylim([0.25,0.45]);
legend('Label 1','Label 2','Label 3','Label 4','Label All');
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_All_zoom_all.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
plot(frames,res_mean1_CEUS(1:L),frames,res_mean2_CEUS(1:L),frames,res_mean3_CEUS(1:L),...
    frames,res_mean4_CEUS(1:L),frames,res_mean_all_CEUS(1:L));
title('TIC (CEUS) - Mean All Labels');
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
%ylim([0.05,0.15]);
legend('Label 1','Label 2','Label 3','Label 4','Label All');
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_All_CEUS_zoom_all.jpg']);
close(f);


% Saves labeled data and results (zoom+normalized):
% B-Mode (zoom+normalized):
f = figure('visible','off');
hold on;
plot(frames,res_mean1(1:L)./res_mean_all(1:L),'LineWidth',3);
errorbar(frames,res_mean1(1:L)./res_mean_all(1:L),res_std1(1:L)./res_std_all(1:L),'LineWidth',1.5);
title(['TIC - Mean + STD, Label 1, ',int2str(sum(img_L(img_L == 1))),...
    ' pixels (~',num2str(pixel_cm^2*sum(img_L(img_L == 1))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_1_zoom_norm.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
plot(frames,res_mean2(1:L)./res_mean_all(1:L),'LineWidth',3);
errorbar(frames,res_mean2(1:L)./res_mean_all(1:L),res_std2(1:L)./res_std_all(1:L),'LineWidth',1.5);
title(['TIC - Mean + STD, Label 2, ',int2str(sum(1/2*img_L(img_L == 2))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/2*img_L(img_L == 2))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_2_zoom_norm.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
plot(frames,res_mean3(1:L)./res_mean_all(1:L),'LineWidth',3);
errorbar(frames,res_mean3(1:L)./res_mean_all(1:L),res_std3(1:L)./res_std_all(1:L),'LineWidth',1.5);
title(['TIC - Mean + STD, Label 3, ',int2str(sum(1/3*img_L(img_L == 3))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/3*img_L(img_L == 3))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_3_zoom_norm.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
plot(frames,res_mean4(1:L)./res_mean_all(1:L),'LineWidth',3);
errorbar(frames,res_mean4(1:L)./res_mean_all(1:L),res_std4(1:L)./res_std_all(1:L),'LineWidth',1.5);
title(['TIC - Mean + STD, Label 4, ',int2str(sum(1/4*img_L(img_L == 4))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/4*img_L(img_L == 4))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_4_zoom_norm.jpg']);
close(f);
%
% CEUS (zoom+normalized):
f = figure('visible','off');
hold on;
plot(frames,res_mean1_CEUS(1:L)./res_mean_all_CEUS(1:L),'LineWidth',3);
errorbar(frames,res_mean1_CEUS(1:L)./res_mean_all_CEUS(1:L),res_std1_CEUS(1:L)./res_std_all_CEUS(1:L),'LineWidth',1.5);
title(['TIC (CEUS) - Mean + STD, Label 1, ',int2str(sum(img_L(img_L == 1))),...
    ' pixels (~',num2str(pixel_cm^2*sum(img_L(img_L == 1))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_1_CEUS_zoom_norm.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
plot(frames,res_mean2_CEUS(1:L)./res_mean_all_CEUS(1:L),'LineWidth',3);
errorbar(frames,res_mean2_CEUS(1:L)./res_mean_all_CEUS(1:L),res_std2_CEUS(1:L)./res_std_all_CEUS(1:L),'LineWidth',1.5);
title(['TIC (CEUS) - Mean + STD, Label 2, ',int2str(sum(1/2*img_L(img_L == 2))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/2*img_L(img_L == 2))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_2_CEUS_zoom_norm.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
plot(frames,res_mean3_CEUS(1:L)./res_mean_all_CEUS(1:L),'LineWidth',3);
errorbar(frames,res_mean3_CEUS(1:L)./res_mean_all_CEUS(1:L),res_std3_CEUS(1:L)./res_std_all_CEUS(1:L),'LineWidth',1.5);
title(['TIC (CEUS) - Mean + STD, Label 3, ',int2str(sum(1/3*img_L(img_L == 3))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/3*img_L(img_L == 3))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_3_CEUS_zoom_norm.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
plot(frames,res_mean4_CEUS(1:L)./res_mean_all_CEUS(1:L),'LineWidth',3);
errorbar(frames,res_mean4_CEUS(1:L)./res_mean_all_CEUS(1:L),res_std4_CEUS(1:L)./res_std_all_CEUS(1:L),'LineWidth',1.5);
title(['TIC (CEUS) - Mean + STD, Label 4, ',int2str(sum(1/4*img_L(img_L == 4))),...
    ' pixels (~',num2str(pixel_cm^2*sum(1/4*img_L(img_L == 4))),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_4_CEUS_zoom_norm.jpg']);
close(f);
%
% All (zoom+normalized):
f = figure('visible','off');
hold on;
plot(frames,res_mean_all(1:L)./res_mean_all(1:L),'LineWidth',3);
errorbar(frames,res_mean_all(1:L)./res_mean_all(1:L),res_std_all(1:L),'LineWidth',1.5);
title(['TIC - Mean + STD, Label All, ',int2str(sum(img_all == img_all)),...
    ' pixels (~',num2str(pixel_cm^2*sum(img_all == img_all)),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_All_zoom_norm.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
plot(frames,res_mean_all_CEUS(1:L),'LineWidth',3);
errorbar(frames,res_mean_all_CEUS(1:L),res_std_all_CEUS(1:L),'LineWidth',1.5);
title(['TIC (CEUS) - Mean + STD, Label All, ',int2str(sum(img_all_CEUS == img_all_CEUS)),...
    ' pixels (~',num2str(pixel_cm^2*sum(img_all_CEUS == img_all_CEUS)),'[cm^2])']);
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_All_CEUS_zoom_norm.jpg']);
close(f);


% Saves figures (images):
% B-Mode:
f = figure('visible','off');
image(uint8(shortvidreg(:,:,:,ind_reg)));
title(['Original B-Mode, frame ',int2str(frame)]);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(gca,xangle);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
xlabel('Lateral [cm]');
ylabel('Axial [cm]');
pause(0.01); % For letting the stupid computer arrange its life.
saveas(f,['.\Results - temp\',description_2,'_bmode.jpg']);
close(f);
%
%{
f = figure('visible','off');
image(img_3);
title(['Segmented B-Mode, frame ',int2str(frame)]);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(gca,xangle);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
xlabel('Lateral [cm]');
ylabel('Axial [cm]');
saveas(f,['.\Results - temp\',description_2,'_bmode_seg.jpg']);
close(f);
%}
%
% CEUS:
f = figure('visible','off');
image(uint8(shortvidreg_2(:,:,:,ind_reg)));
title(['Original CEUS, frame ',int2str(frame)]);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(gca,xangle);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
xlabel('Lateral [cm]');
ylabel('Axial [cm]');
pause(0.01); % For letting the stupid computer arrange its life.
saveas(f,['.\Results - temp\',description_2,'_ceus.jpg']);
close(f);
%
%{
f = figure('visible','off');
image(img_3_CEUS);
title(['Segmented CEUS, frame ',int2str(frame)]);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(gca,xangle);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
xlabel('Lateral [cm]');
ylabel('Axial [cm]');
saveas(f,['.\Results - temp\',description_2,'_ceus_seg.jpg']);
close(f);
%}


% Plots with smoothing:
f = figure('visible','off');
hold on;
plot(frames,res_mean1_pol(1:L),frames,res_mean2_pol(1:L),frames,res_mean3_pol(1:L),...
    frames,res_mean4_pol(1:L),frames,res_mean_all_pol(1:L));
title('TIC - Mean All Labels smoothed');
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
%ylim([0.25,0.45]);
legend(['Label 1 - Avg ',num2str(mean(res_mean1)),', STD ',num2str(mean(res_std1))],...
    ['Label 2 - Avg ',num2str(mean(res_mean2)),', STD ',num2str(mean(res_std2))],...
    ['Label 3 - Avg ',num2str(mean(res_mean3)),', STD ',num2str(mean(res_std3))],...
    ['Label 4 - Avg ',num2str(mean(res_mean4)),', STD ',num2str(mean(res_std4))],...
    ['Label All - Avg ',num2str(mean(res_mean_all)),', STD ',num2str(mean(res_std_all))]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_All_zoom_all_avg.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
plot(frames,res_mean1_CEUS_pol(1:L),frames,res_mean2_CEUS_pol(1:L),frames,res_mean3_CEUS_pol(1:L),...
    frames,res_mean4_CEUS_pol(1:L),frames,res_mean_all_CEUS_pol(1:L));
title('TIC (CEUS) - Mean All Labels smoothed');
xlabel('Frames');
ylabel('Intensity ([0,1])');
xlim([0,L+1]);
%ylim([0.05,0.15]);
legend(['Label 1 - Avg ',num2str(mean(res_mean1_CEUS)),', STD ',num2str(mean(res_std1_CEUS))],...
    ['Label 2 - Avg ',num2str(mean(res_mean2_CEUS)),', STD ',num2str(mean(res_std2_CEUS))],...
    ['Label 3 - Avg ',num2str(mean(res_mean3_CEUS)),', STD ',num2str(mean(res_std3_CEUS))],...
    ['Label 4 - Avg ',num2str(mean(res_mean4_CEUS)),', STD ',num2str(mean(res_std4_CEUS))],...
    ['Label All - Avg ',num2str(mean(res_mean_all_CEUS)),', STD ',num2str(mean(res_std_all_CEUS))]);
set(gcf,'Position',get(0,'Screensize'));
hold off;
saveas(f,['.\Results - temp\',description_2,'_TIC_All_CEUS_zoom_all_avg.jpg']);
close(f);


% Histograms plot:
% B-Mode:
f = figure('visible','off');
hold on;
h = histogram(res_hist1,(0:1:256)-0.5,'Normalization','probability');
plot(0:1:255,h.Values)
hold off;
title('histogram of Label 1');
saveas(f,['.\Results - temp\',description_2,'_hist1.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
h = histogram(res_hist2,(0:1:256)-0.5,'Normalization','probability');
plot(0:1:255,h.Values)
hold off;
title('histogram of Label 2');
saveas(f,['.\Results - temp\',description_2,'_hist2.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
h = histogram(res_hist3,(0:1:256)-0.5,'Normalization','probability');
plot(0:1:255,h.Values)
hold off;
title('histogram of Label 3');
saveas(f,['.\Results - temp\',description_2,'_hist3.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
h = histogram(res_hist4,(0:1:256)-0.5,'Normalization','probability');
plot(0:1:255,h.Values)
hold off;
title('histogram of Label 4');
saveas(f,['.\Results - temp\',description_2,'_hist4.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
h = histogram(res_hist_all,(0:1:256)-0.5,'Normalization','probability');
plot(0:1:255,h.Values)
hold off;
title('histogram of Label All');
saveas(f,['.\Results - temp\',description_2,'_hist_all.jpg']);
close(f);
%
% CEUS:
f = figure('visible','off');
hold on;
h = histogram(res_hist1_CEUS,(0:1:256)-0.5,'Normalization','probability');
plot(0:1:255,h.Values)
hold off;
title('histogram of Label 1 (CEUS)');
saveas(f,['.\Results - temp\',description_2,'_hist1_CEUS.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
h = histogram(res_hist2_CEUS,(0:1:256)-0.5,'Normalization','probability');
plot(0:1:255,h.Values)
hold off;
title('histogram of Label 2 (CEUS)');
saveas(f,['.\Results - temp\',description_2,'_hist2_CEUS.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
h = histogram(res_hist3_CEUS,(0:1:256)-0.5,'Normalization','probability');
plot(0:1:255,h.Values)
hold off;
title('histogram of Label 3 (CEUS)');
saveas(f,['.\Results - temp\',description_2,'_hist3_CEUS.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
h = histogram(res_hist4_CEUS,(0:1:256)-0.5,'Normalization','probability');
plot(0:1:255,h.Values)
hold off;
title('histogram of Label 4 (CEUS)');
saveas(f,['.\Results - temp\',description_2,'_hist4_CEUS.jpg']);
close(f);
%
f = figure('visible','off');
hold on;
h = histogram(res_hist_all_CEUS,(0:1:256)-0.5,'Normalization','probability');
plot(0:1:255,h.Values)
hold off;
title('histogram of Label All (CEUS)');
saveas(f,['.\Results - temp\',description_2,'_hist_all_CEUS.jpg']);
close(f);
end