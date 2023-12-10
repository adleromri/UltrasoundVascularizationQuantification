%% Rings plot
%
% Syntax:
% PanGUI_RingsPlot(img_3,img_3_CEUS,frame,xangle,xticks,xticklabels,yticks,yticklabels,respar)
%
% Input:
% ...

function PanGUI_RingsPlot(img_3,img_3_CEUS,frame,xangle,xticks,xticklabels,yticks,yticklabels,respar,description_0)
% B-Mode:
f = figure('visible','off');
image(img_3);
title(['Segmented B-Mode, frame ',int2str(frame)]);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(gca,xangle);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
xlabel('Lateral [cm]');
ylabel('Axial [cm]');
saveas(f,[respar,'bmode_seg',description_0,'.jpg']);
close(f);

% CEUS:
f = figure('visible','off');
image(img_3_CEUS);
title(['Segmented CEUS, frame ',int2str(frame)]);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(gca,xangle);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
xlabel('Lateral [cm]');
ylabel('Axial [cm]');
saveas(f,[respar,'ceus_seg',description_0,'.jpg']);
close(f);
end