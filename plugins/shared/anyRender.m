function [srim,hax]=anyRender(locD,p,varargin)
% no arguments: returns required parameters.
%locD,p
%removefilter (fields)
%'x','y', 'xerr','yerr'

%also determine beforehand if grouped or not

if nargin==0
    srim={'layers','sr_pixrec','sr_pos','sr_size','numberOfLayers','sr_layerson','cam_pixelsize_nm','sr_roihandle','linewidth_roi'};
    return;
end

pp=roiparser(varargin);

if ~isfield(p,'sr_axis')
    p.sr_axis=[];
end

if isempty(pp.rangex)
    pp=rmfield(pp,'rangex');
end
if isempty(pp.rangey)
    pp=rmfield(pp,'rangey');
end

if isempty(pp.xpix)
    pp.xpix=p.sr_pixrec;
end
if isempty(pp.ypix)
    pp.ypix=p.sr_pixrec;
end
p.sr_pixrec=[pp.xpix,pp.ypix];

 
% lochere=locD.copy();
if ~isempty(pp.within)
    indroi=pp.within;
else
[~,indroi]=locD.getloc('xnm','position',pp.position);

end


% fields={pp.x,pp.y,pp.sx,pp.sy,'xnm','ynm','znm','intensity_render'};
% for k=1:p.numberOfLayers
%     fields{end+1}=p.(['layer' num2str(k) '_']).renderfield.selection;
% end
% clearfilter=false;
% lochere=locD.copy(fields,indroi,clearfilter);
lochere=locD.copy;


%later: if sx or sy scalar: use for all
if strcmp(pp.x,'xnmline')||strcmp(pp.y,'xnmline')||strcmp(pp.x,'ynmline')||strcmp(pp.y,'ynmline')
   if pp.groupstate(1)
    [locl,indin]=lochere.getloc({'xnmline','ynmline'},'position','roi','grouping','ungrouped');
    locline=zeros(length(indin),1);
    locline(indin)=locl.xnmline;
    lochere.loc.xnmline=locline;
    locline=zeros(length(indin),1);
    locline(indin)=locl.ynmline;
    lochere.loc.ynmline=locline;
   end
   if pp.groupstate(2)
    [locl,indin]=lochere.getloc({'xnmline','ynmline'},'position','roi','grouping','grouped');
    locline=zeros(length(indin),1);
    locline(indin)=locl.xnmline;
    lochere.grouploc.xnmline=locline;
    locline=zeros(length(indin),1);
    locline(indin)=locl.ynmline;
    lochere.grouploc.ynmline=locline;
   end
    
    roih=p.sr_roihandle;
    if isvalid(roih)
        rpos=roih.getPosition;
    %     mpos=mean(rpos,1);
        lr=sqrt(sum((rpos(2,:)-rpos(1,:)).^2));
        rx=[-lr/2 lr/2]*1000;
        ry=[-0.5 0.5]*p.linewidth_roi;
    else
    rx=[min(locs.xnmline) max(locs.xnmline)];
    ry=[min(locs.ynmline) max(locs.ynmline)];
    end
    if strcmp(pp.x,'xnmline')
        pp.rangex=rx;
    end
    
    if strcmp(pp.x,'ynmline')
        pp.rangex=ry;
    end
    if strcmp(pp.y,'xnmline')
        pp.rangey=rx;
    end
    
    if strcmp(pp.y,'ynmline')
        pp.rangey=ry;
    end
        
end

if pp.groupstate(1)
lochere.loc.x=lochere.loc.(pp.x);
lochere.loc.y=lochere.loc.(pp.y);
v1=lochere.loc.(pp.sx);
v2=lochere.loc.(pp.sy);
lochere.loc.sx=v1;
lochere.loc.sy=v2;
end
% lochere.loc.xnm=v1;
% lochere.loc.ynm=v2;

if pp.groupstate(2)
lochere.grouploc.x=lochere.grouploc.(pp.x);
lochere.grouploc.y=lochere.grouploc.(pp.y);
v1=lochere.grouploc.(pp.sx);
v2=lochere.grouploc.(pp.sy);
lochere.grouploc.sx=v1;
lochere.grouploc.sy=v2;
end
% renderer.setParameters(p,pp);
for k=1:p.numberOfLayers
    pl=p.(['layer' num2str(k) '_']);
    
    if pl.layercheck
        if isempty(strfind(pp.x,'nm'))||isempty(strfind(pp.y,'nm'))
       
            pl.mingaussnm=0;%why????
        end
         pr=copyfields(copyfields(pl,p),pp);
         if pl.groupcheck
             indroi=lochere.getloc('ingrouped','layer',k,'position',pp.position).ingrouped;  
            lochere.layer(k).images.srimage=renderSMAP(lochere.grouploc,pr,k,indroi);
         else
             indroi=lochere.getloc('inungrouped','layer',k,'position',pp.position).inungrouped;  
             lochere.layer(k).images.srimage=renderSMAP(lochere.loc,pr,k,indroi);
         end

        lochere.layer(k).images.srimage.rangex=lochere.layer(k).images.srimage.rangex;
        lochere.layer(k).images.srimage.rangey=lochere.layer(k).images.srimage.rangey;
        lochere.layer(k).images.finalImages=drawerSMAP(lochere.layer(k).images.srimage,pr);        
 
    end
end
srim=displayerSMAP(lochere.layer,pr);

end

function pres=roiparser(args)
% fields{end+1}='all';
p = inputParser;   
% addOptional(p,'fields',{},@(x)  any(myvalidatestring(x,fields)));
addParameter(p,'removefilter',{});
addParameter(p,'x','xnm');
addParameter(p,'y','ynm');
addParameter(p,'sx','locprecnm');
addParameter(p,'sy','locprecnm');
addParameter(p,'xpix',[]);
addParameter(p,'ypix',[]);
addParameter(p,'rangex',[]);
addParameter(p,'rangey',[]);
addParameter(p,'position','all');
addParameter(p,'within',[]);
addParameter(p,'groupstate',[1 1]);
parse(p,args{:});
pres=p.Results;
end
