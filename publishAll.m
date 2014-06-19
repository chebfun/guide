%%
% NB that you will have to manuall put a "width=6in" into guide10.tex for
% the image of the chebgui; otherwise it is much too large.


%%
guideFiles = dir('guide*.m');
opt.format = 'latex';
opt.outputDir = '/Users/driscoll/Documents/guideLatex';
opt.imageFormat = 'epsc'; 
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
opt.outputDir = '/Users/driscoll/Documents/guideLatex';
opt.stylesheet = 'publish2latex.xsl';
for j = 1:length(guideFiles)
    defaultSettings
    publish(guideFiles(j).name,opt);
end
