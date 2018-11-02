function pnew=convertparameters(pold)

replace={'string','String';'value','Value';'sitefov','se_sitefov';'cellfov','se_cellfov';...
    'sitepixelsize','se_sitepixelsize';'cellpixelsize','se_cellpixelsize';'siteroi','se_siteroi';...
    'rotate','se_rotate';'drawboxes','se_drawboxes'};

pnew=removepar(pold);

s=size(replace);
for k=1:s(1)
pnew=replacefieldsrecursive(pnew,replace{k,1},replace{k,2});
end


function pnew=removepar(pold)
pnew=[];
if isfield(pold,'par')
    pnew=copyfields(pnew,pold.par);
end
if isfield(pold,'children')
    fn=fieldnames(pold.children);
    for k=1:length(fn);
    pnew.children.(fn{k})=removepar(pold.children.(fn{k}));
    end
end

