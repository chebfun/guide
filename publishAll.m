%%
% NB that you will have to manually put a "width=6in" into guide10.tex for
% the image of the chebgui; otherwise it is much too large.
addpath bookstuff

%%
guideFiles = dir('guide*.m');
opt.format = 'latex';
opt.outputDir = 'latex/color';
opt.imageFormat = ''; 
opt.stylesheet = 'publish2latex.xsl';
for j = 1:length(guideFiles)
    defaultSettings
    publish(guideFiles(j).name,opt);
end

%%
% Black and white
guideFiles = dir('guide*.m');
opt.format = 'latex';
opt.imageFormat = 'eps'; 
opt.outputDir = 'latex/bw';
opt.stylesheet = 'publish2latex.xsl';
for j = 1:length(guideFiles)
    defaultSettings
    publish(guideFiles(j).name,opt);
end
