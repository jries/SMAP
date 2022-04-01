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

        %%%IC220401 - assign all model params based on the 3D polyline
        if strcmp(convertCharsToStrings(target),'usr_all') && contains(rule,'{ALLTARGETS}')
                        
            initpars=eval(erase(rule,'{ALLTARGETS}')); %get source cell 
            ruleExprs=rule; %save common rule
            matchStr = regexp(ruleExprs,'obj\.converterSource\{\d+\}','match'); cs=str2num(matchStr{1}(end-1)); %get converterSource number
            matchStr = regexp(ruleExprs,'(usr|post)_\w+','match'); parsVar=matchStr{1}; %get the name of the parameter containing variable (not tested with post_)


            for p=1:size(initpars,2)
                rule=regexprep(ruleExprs, '(\{ALLTARGETS\})', ['\)\.value\(strcmp\(\{obj\.converterSource\{' num2str(cs) '\}\.getVariable\(' char(39) parsVar char(39) '\)\.name\}, ' char(39) replace(initpars(p).name,'.','\.') char(39) '\)\)']);
                rule=['struct2table(' rule];
                %mulitpletargets=initpars(p).name;
                target=initpars(p).name;
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
            
        else % insert original code (value=...end of IC220401 addition)
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
            end %%% end of IC220401 addition
        end 
        
        
    end
end
end