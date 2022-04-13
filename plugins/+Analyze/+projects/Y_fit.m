function out=Y_fit(varargin)
newpath=strrep(pwd,'SMAP','ries-private');
if exist(newpath,'dir')
    if ~isdeployed
        addpath([newpath filesep 'VectorPSF_Fit']);
    end
    out=Y_fit_plugin(varargin{:});
else
    out=[];
end