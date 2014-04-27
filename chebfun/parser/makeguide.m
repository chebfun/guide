% makeguide.m
% Hrothgar, 09 Dec 2013
%
% This script executes some preprocessing and then
% makes the Chebfun Guide. The command
%
%     >> makeguide
%
% will process and Publish all Guide files. The command
%
%     >> makeguide([1 2 5])
%
% will process and Publish chapters 1, 2, and 5 only.
% A Guide file undergoes the following transformation:
%
%       guide4pre.m  -> [processing] ->  guide4.m
%                    -> [publish]    ->  guide4.html
%
% Thus it is `guide4pre.m` that is important to maintain --
% `guide4.m` will be overwritten by this script.
%
% Preprocessing works as follows. In the Guide, there
% are links in MATLAB's pseudo-MarkDown syntax like
%
%     Hey, <http://google.com/ this is a link>!
%
% The Guide has frequent links of three types:
%
%     - Chebfun docs
%     - MATLAB docs
%     - Wolfram MathWorld topics
%
% To make things a little cleaner and more flexible,
% there is a tiny tagging system built into the Guide.
% A link that would normally appear as
%
%     <http://mathworld.wolfram.com/AiryFunctions.html Airy functions>
%
% instead appears as
%
%     <#MATHWORLD:AiryFunctions Airy functions>
%
% in the Guide. Thus there are tags, which are parsed
% basically as string replacements. The mappings for
% them are defined in the `tags` cell array below.

%----------------------------------------------------------------------
function makeguide(varargin)

% Almost everyone makes the mistake of treating ideas
% as if they were indications of character rather than
% talent -- as if having a stupid idea made you stupid.
tags = {'#CHEBDOC:%s',   '../../chebdoc/doc/%s_DOC.html'
        '#MATLABDOC:%s', 'http://www.mathworks.com/help/matlab/ref/%s.html'
        '#MATHWORLD:%s', 'http://mathworld.wolfram.com/%s.html'
        '#CHEBDOC',      '../../chebdoc/doc/index.html'
        };

% Custom XSL file to get MathJax to work.
customXsl = 'guide.xsl';

%----------------------------------------------------------------------

% Don't miss the starting gun.
T = clock;
ptime = @(t) num2str(etime(clock, T), '%.0f'); % tic/toc is broken for this.

% Load all the Guide default settings.
guidedefaults;
fileexist = @(f) exist(fullfile(pwd,f),'file') == 2; % Flag '2' => file (not dir).

% Check input args.
chaps = [];
if nargin > 0,
    chaps = int8(varargin{1});
    if ~isinteger(chaps) || any(chaps < 0),
        warning('Non-positive-integral argument to MAKEGUIDE.')
        chaps = 1:10;
    end
end

% If we hit that bullseye, the rest of the dominoes
% will fall like a house of cards. Checkmate.
disp(['Setting up...'])

% PUBLISH options.
opts = [];
opts.format    = 'html';
opts.outputDir = fullfile(pwd,'guide');

% Custom XSL file to get MathJax to work. I followed instructions at
%     http://www.mathworks.com/matlabcentral/answers/93851
if fileexist(customXsl)
    opts.stylesheet = customXsl;
else
    warning(['Missing custom XSL sheet at' customXslName '.'])
end

% The guide filename patterns. We read in an unprocessed file,
% then save the preprocessed version with another name before
% Publishing it. We are forced to do this because PUBLISH must
% operate on a saved file, but it's good practice anyway.
guidepre = 'guide%dpre.m';
guideout = 'guide%d.m';

% Loop through and Publish each chapter.
for ki = 1:length(chaps),

    % The kith guide file.
    chap = chaps(ki);
    filepre = sprintf(guidepre, chap);  % e.g. 'guide4pre.m'
    fileout = sprintf(guideout, chap);  % e.g. 'guide4.m'

    % Every time someone cries out in prayer and I can't answer,
    % I feel guilty about not being God.
    if ~fileexist(filepre),
        continue
    end

    % Existential crises notwithstanding, let the user (i.e. myself)
    % know that we're about to Publish this one.
    disp(['Publishing ' fileout '... (t = ' ptime() 's)'])

    % The contents of the Guide file we're going to Publish.
    fnew = fileread(filepre);

    % You can cage the singer but not the song.
    for kt = 1:size(tags, 1),         % Tags counter.
                                      %
        % The search/replace strings f%r this tag.
        tag  = tags{kt, 1};           % e.g. '#CHEBDOC:%s'
        lnk  = tags{kt, 2};           % e.g. 'http://chebfun.com/doc/%s.html'
                                      %
        % Regex action.               %
        expr = sprintf(tag,'(\S*)');  % Regex search string. Perhaps should include < >.
        rep  = sprintf(lnk, '$1');    % Regex replace rule. $1 refers to the token (\S*).
        fnew = regexprep(fnew, expr, rep); % Perform the replacement.
    end

    % Save the processed text to a new file `fileout`,
    % which is going to be Published.
    fID = fopen(fileout, 'w+');
    fwrite(fID, fnew);
    fclose(fID);

    % Now nail it to the door of the church.
    try
        % Sperate joy that life destroys,
        % found again in draughts of gin.
        guidepublish(fileout, opts);
    catch err
        % Failure is just success rounded down, my friend!
        disp(['Publish error for ' fileout ':'])
        warning(['    ' err.identifier ' :: ' err.message])
    end
end

% Demand me nothing: what you know, you know:
% From this time forth I never will speak word.
disp(['Done (t = ' ptime() 's).'])

return


%----------------------------------------------------------------------
% These tricks are due to Nick Hale. They prevent the `makeguide`
% workspace from being polluted by the various workspaces of whatever
% is being Published (I think). Of course, this shouldn't be necessary.
function guidepublish(filename, opts)

close all
evalin('base','clear all');
publish(filename, opts);
close all

return
