
%%
guideFiles = dir('guide*.m');
opt.format = 'latex';
opt.outputDir = '../guideLatex';
%opt.figureSnapMethod = 'print';
opt.stylesheet = 'publish2latex.xsl';
for j = 1:length(guideFiles)
    defaultSettings
    publish(guideFiles(j).name,opt);
end
