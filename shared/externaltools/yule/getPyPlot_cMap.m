function cmp = getPyPlot_cMap(nam,n,keepAlpha,pyCmd)
% cmp = getPyPlot_cMap(nam [, n, keepAlpha, pyCmd])
%
%
% ::INPUT::
%
% nam:      Colormap name from matplotlib library (case sensitve!). See
%           below. Alternatively, '!GetNames' returns a cellstring of
%           available colormaps.
% n:        Number of Colors; defaults to 128
% keepAlpha: Switch to keep the alpha channel of the colormap (4th colum);
%           defaults to false. If true, a Nx4 matrix is returned in cmp
%           (instead of Nx3).
% pyCmd:    python command; defaults to 'python'
%
% 
% ::OUTPUT::
%
% cmp       A Nx3 (Nx4 if keepAlpha is true) double array of RGB(A) values.
%
%
% Colormap name can be one of the following:
%
%     Accent      gist_earth        Oranges           RdYlBu
%     Accent_r    gist_earth_r      Oranges_r         RdYlBu_r
%     afmhot      gist_gray         OrRd              RdYlGn
%     afmhot_r    gist_gray_r       OrRd_r            RdYlGn_r
%     autumn      gist_heat         Paired            Reds
%     autumn_r    gist_heat_r       Paired_r          Reds_r
%     binary      gist_ncar         Pastel1           seismic
%     binary_r    gist_ncar_r       Pastel1_r         seismic_r
%     Blues       gist_rainbow      Pastel2           Set1
%     Blues_r     gist_rainbow_r    Pastel2_r         Set1_r
%     bone        gist_stern        pink              Set2
%     bone_r      gist_stern_r      pink_r            Set2_r
%     BrBG        gist_yarg         PiYG              Set3
%     BrBG_r      gist_yarg_r       PiYG_r            Set3_r
%     brg         GnBu              PRGn              Spectral
%     brg_r       GnBu_r            PRGn_r            spectral
%     BuGn        gnuplot           prism             spectral_r
%     BuGn_r      gnuplot_r         prism_r           Spectral_r
%     BuPu        gnuplot2          PuBu              spring
%     BuPu_r      gnuplot2_r        PuBu_r            spring_r
%     bwr         gray              PuBuGn            summer
%     bwr_r       gray_r            PuBuGn_r          summer_r
%     CMRmap      Greens            PuOr              terrain
%     CMRmap_r    Greens_r          PuOr_r            terrain_r
%     cool        Greys             PuRd              winter
%     cool_r      Greys_r           PuRd_r            winter_r
%     coolwarm    hot               Purples           Wistia
%     coolwarm_r  hot_r             Purples_r         Wistia_r
%     copper      hsv               rainbow           YlGn
%     copper_r    hsv_r             rainbow_r         YlGn_r
%     cubehelix   jet               RdBu              YlGnBu
%     cubehelix_r jet_r             RdBu_r            YlGnBu_r
%     Dark2       nipy_spectral     RdGy              YlOrBr
%     Dark2_r     nipy_spectral_r   RdGy_r            YlOrBr_r
%     flag        ocean             RdPu              YlOrRd
%     flag_r      ocean_r           RdPu_r            YlOrRd_r
% 
% V 1.3; Konrad Schumacher, 07.2019

if strcmpi(nam,'!GetNames')
    % translate switch to retrieve colormap names into python-switch:
    nam = 'listCMapNames';
end


% defaults:
if ~exist('n','var') || isempty(n)
    n = 128;
end
if ~exist('keepAlpha','var') || isempty(keepAlpha)
    keepAlpha = 0;
end
if ~exist('pyCmd','var') || isempty(pyCmd)
    pyCmd = 'python';
end


% check if python script is present
pyScript = which('pyplotCMap2txt.py');
assert(~isempty(pyScript), 'getPyPlot_cMap:PyScriptNotFound', ...
    'Could not find python script (%s).','pyplotCMap2txt.py');

tmpf = tempname;

% call python script
comd = sprintf('%s "%s" %s -o "%s" -n %d', pyCmd, pyScript, nam, tmpf, n);
[s,m] = system(comd);

% check if system command ran w/o error
assert(s==0, 'getPyPlot_cMap:SystemCMDFailed', ...
        'There was an error executing the command\n\t%s\nSystem returned:\n\t%s', ...
        comd, m);


if strcmp(nam,'listCMapNames')
    % cMap names retrieved; read text file
    fid = fopen(tmpf,'r');
    cmp = textscan(fid,'%s');
    fclose(fid);
    cmp = cmp{1};
    
else
    % load cMap data from text file
    cmp = load(tmpf,'-ascii');

    if keepAlpha
    else % remove 4th column of alpha values
        cmp = cmp(:,1:3);
    end
end


% delete temp-file
delete(tmpf);

end

