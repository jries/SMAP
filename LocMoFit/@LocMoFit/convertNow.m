function convertNow(obj,locs,varargin)
%% CONVERTNOW Evaluate parameter conversion right away
% Conver locs' values or previously fitted parameters to the initial guesses
% or temporary variables.
% Skip if there is no any converterRules.

%% Deal with the varargin
p = inputParser;
p.addParameter('locs2', [])
p.parse(varargin{:})
locs2 = p.Results.locs2;                % the second set of locs

%% Main
% Skip when converterRules is empty:
if ~isempty(obj.converterRules)
    for k = 1:length(obj.converterRules.rule)
        target = obj.converterRules.target{k};
        rule = obj.converterRules.rule{k};
        value = eval(rule);
        if startsWith(target,'usr_')||startsWith(target,'post_')
            obj.converterUserDefined.(target)=value;
        else
            % check it is which of the following:
            %  1) m1.lPar.[fn]
            %  2) m1.lPar.[fn].[lb]
            targetParts = strsplit(target,'.');
            switch length(targetParts)
                
                case 3
                    % for the case (1) above
                    obj.setParArg(target,'value',value)
                case 4
                    % for the case (2) above
                    obj.setParArg(target,targetParts{end},value)
            end
        end
    end
end
end