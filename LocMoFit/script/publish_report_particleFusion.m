%% Publish reports of particle fusion analyses
% One can use this script to generate reports of their particle fusion
% analyses

%% Specify the folder where the results of template free particle fusion are
rootPath = 'D:\path\to\the\folder\where\the\output\files\are\saved\';

% Example:
% rootPath = 'D:\particleFusion\analysis\';

%% ---------Typically, users do not have to change the following code------------
%% Prepare for generating the report

[filepath,name] = fileparts(mfilename('fullpath'));
filepath = [filepath filesep];

%% Stop recording the messages in the MATLAB command window
set(0,'DefaultFigureVisible','off')
echo off;
warning off

options = [];
options.showCode = false;
options.format = 'html';
options.maxOutputLines = 0;

%% Detect all the analyses and generate reports for them
fileList = dir([rootPath '*avg.mat']);

for k = 1:length(fileList)
    fName = fileList(k).name;
    ID_analysis = regexprep(fName,'_avg.mat$','');
    publish([filepath 'generate_report_particleFusion.m'],options)
    if exist([filepath ID_analysis], 'dir')
        rmdir([filepath ID_analysis], 's');
    end
    movefile([filepath 'html'], [filepath ID_analysis])
end

%% Reset the message display
pause(5)
set(0,'DefaultFigureVisible','on')
echo on;
warning on;
disp('Done');