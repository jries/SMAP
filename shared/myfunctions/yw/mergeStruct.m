function s = mergeStruct(varargin)
    % it requires the input struct to have the same fields
    fn = fieldnames(varargin{1});
    s = [];
    toRun1 = cellstr(append('varargin{', string(num2str((1:nargin)')),'}.(fn{k}); '));
    toRun2 = cellstr(append('varargin{', string(num2str((1:nargin)')),'}.(fn{k}) '));
    for k = 1:length(fn)
        try
            s.(fn{k}) = eval(['[' toRun1{:} ']']);
        catch
            s.(fn{k}) = eval(['[' toRun2{:} ']']);
        end
    end
end