function optionFeedback = createOptionFeedback(options)
%CREATEOPTIONFEEDBACK check if options are default and create structure used
%by exit messages for Optimization Toolbox solvers.
%
% This utility takes the input options structure and makes a "mirror"
% structure (optionFeedback) that contains fields for the options used in
% the exit messages for Optimization Toolbox solvers. 
%
% The fields of optionFeedback contain the strings 'default' or
% 'selected' depending on whether or not the user set the
% corresponding option. The strings in the structure are then used by
% createExitMsg to give feedback to the user on whether the options that
% affect the stopping condition had been changed.

%   Copyright 2008 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2008/12/01 07:40:32 $

% Make sure that options contains all fields
options = optimset(options);

% Check if TolFun is the default value
if isempty(options.TolFun)
    optionFeedback.TolFun = 'default';
else
    optionFeedback.TolFun = 'selected';
end

% Check if TolX is the default value
if isempty(options.TolX)
    optionFeedback.TolX = 'default';
else
    optionFeedback.TolX = 'selected';
end

% Check if MaxIter is the default value
if isempty(options.MaxIter)
    optionFeedback.MaxIter = 'default';
else
    optionFeedback.MaxIter = 'selected';
end

% Check if MaxFunEvals is the default value
if isempty(options.MaxFunEvals)
    optionFeedback.MaxFunEvals = 'default';
else
    optionFeedback.MaxFunEvals = 'selected';
end

% Check if TolCon is the default value
if isempty(options.TolCon)
    optionFeedback.TolCon = 'default';
else
    optionFeedback.TolCon = 'selected';
end

% Check if MaxSQPIter is the default value
if isempty(options.MaxSQPIter)
    optionFeedback.MaxSQPIter = 'default';
else
    optionFeedback.MaxSQPIter = 'selected';
end

% Check if ObjectiveLimit is the default value
if isempty(options.ObjectiveLimit)
    optionFeedback.ObjectiveLimit = 'default';
else
    optionFeedback.ObjectiveLimit = 'selected';
end
