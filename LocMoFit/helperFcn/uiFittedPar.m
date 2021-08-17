function hGui_par = uiFittedPar(varargin)
%
% Usuage:
%    uiFittedPar(ax, fitter)
%    uiFittedPar(fig, fitter)
% Parameters:
%    ax: a MATLAB axes object.
%    fig: a MATLAB figure object.
%    fitter: a LocMoFit object.

if isa(varargin{1},'matlab.graphics.axis.Axes')
    target = varargin{1}.Parent;
else
    target = varargin{1};
end
fitter = varargin{2};
fittedPar = findobj(target, 'Style', 'edit');
FP = fitter.viewPars;
if isempty(fittedPar)
    hGui_par.fittedPar = uitable(target,...
        'Position',[1 1 5.2 11]);
    hGui_par.copy_FP = uicontrol(target,...
        'Style', 'pushbutton',...
        'String', 'Copy all',...
        'Position',[1 12 1 1],...
        'Callback', {@copyTable_callback, FP});
    guiStyle(hGui_par, {'fittedPar', 'copy_FP'})
else
    hGui_par.fittedPar = fittedPar;
end

header = FP.Properties.VariableNames;
hGui_par.fittedPar.ColumnName = header;
hGui_par.fittedPar.Data = table2cell(FP);
end