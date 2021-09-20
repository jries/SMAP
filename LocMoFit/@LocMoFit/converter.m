function converter(obj, source, rule, target)
% converter saves the rule for converting info to the target
% parameter.
%% deal with source
if isempty(obj.converterSource)
    % initiate the converterSource
    obj.converterSource{1}=obj;             % The obj itself will always be at the first position
    k = 1;
    if ~isempty(source)
        obj.converterSource{2}=source;
        k = 2;
    end
else
    if isempty(source)
        % If the source is empty, refer the source to the obj itself, which is at the first position
        k = 1;
    else
        % Otherwise check whether the source is there or not.
        % If it is not there, then add it after the last
        % position.
        l = 1;
        while l <= length(obj.converterSource)
            if isequal(source,obj.converterSource{l})
                k = l;
                break
            else
                k = 0;
            end
            l = l+1;
        end
        isequal({source},obj.converterSource);       % the kth position in the obj.converterSource
        if k == 0
            obj.converterSource{end+1} = source;
        else
        end
    end
end

%% complete the rule
if isempty(rule)
    return
end
ruleExprs = regexprep(rule, '(pars\.m\d+\.[ml]Par\.\w+)', ['obj\.converterSource\{' num2str(k) char("\}\.getVariable(\'$1\')")]);          % Parameters
ruleExprs = regexprep(ruleExprs, '(pars\.m\d+\.offset\.\w+)', ['obj\.converterSource\{' num2str(k) char("\}\.getVariable(\'$1\')")]);      % Offset
ruleExprs = regexprep(ruleExprs, '(usr\_\w+)\((\d+)\)', 'vectorSubset\($1\,$2\)');                                                      % User defined (vector)
ruleExprs = regexprep(ruleExprs, '(usr\_\w+)', ['obj\.converterSource\{' num2str(k) char("\}\.getVariable(\'$1\')")]);                     % User defined (single element)
ruleExprs = regexprep(ruleExprs, '(rel\(\w+)', ['obj\.converterSource\{' num2str(k) char("\}\.$1")]);                                   % Relative info
ruleExprs = regexprep(ruleExprs, 'find\.(m\d+\.\w+)', ['obj\.converterSource\{' num2str(k) char("\}\.getVariable(\'$1\')")]);           % Search for the values
ruleExprs = regexprep(ruleExprs, 'find\.(\w+)', ['obj\.converterSource\{' num2str(k) char("\}\.getVariable(\'$1\')")]);           % Search for the values
%% deal with the rule
if isempty(obj.converterRules)
    % initiate the converterRules
    obj.converterRules.target = {};
    obj.converterRules.rule = {};
end

%% Check the rule for the target is set or not.
% One target should always have only one rule.

if startsWith(target, 'post_')
    if ~isfield(obj.converterUserDefined, target)
        % initiate the converterUserDefined
        obj.converterUserDefined.(target) = 0;
    else
        allStartWith = startsWith(fieldnames(obj.converterUserDefined), target);
        target = [target '_' num2str(sum(allStartWith)+1)];
        obj.converterUserDefined.(target) = 0;
        obj.converterRules.target{end+1} = target;
        obj.converterRules.rule{end+1} = ruleExprs;
    end
end

[~,ind] = ismember(target,obj.converterRules.target);
if ~(ind == 0)
    % If the target is set, then overwrite it
    obj.converterRules.target{ind} = target;
    obj.converterRules.rule{ind} = ruleExprs;
else
    % Otherwise add the rule after the last position.
    obj.converterRules.target{end+1} = target;
    obj.converterRules.rule{end+1} = ruleExprs;
end

if startsWith(target, 'usr_')
    target = replace(target,'usr_','');
    if ~isfield(obj.converterUserDefined, target)
        % initiate the converterUserDefined
        obj.converterUserDefined.(target) = 0;
    end
end
end