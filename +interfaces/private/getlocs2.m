function [locsout,indcombined,hroio]=getlocs2(locData,fields,varargin)
%now each layer is treated separately

%returns subset of localizations. 
%[locs,handle to roi]=getloc(obj,{fields},'PropertyName','Property')
%fields: any locData field, in addition: 'ingrouped',
%'inungrouped','within'
% inlayeru, inlayerg: cell array with indices for each layer

%Properties:
%'grouping': 'grouped', 'ungrouped' (second default),'layer'(default). If not set and combined with
%'layer, use individual layer grouping. If any of layers is ungrouped,
%return ungrouped localizations.
%'channels': double vector of channels
%'filenumber': double vector of filenumbers
%'position': 'all' (default),'roi','fov': position filter
%'position': vector: [centerx, centery, widhtx widthy] or [centerx,
%centery,radius] for circular ROI
%   'interfaces.SEsites: locaizations in site.
%'layer': double number or vector of layers.

%'removeFilter': cell array of filter names to remove 
%'within', indices which localizations to consider.
%'shiftxy', shifts x, y  and z by specified vector
%'locselector',p structure from gui


%no localizations present, return empty
% indlayer={};
if isempty(locData.loc)||isempty(locData.loc.(getelement(fieldnames(locData.loc),1)))
    for k=1:length(p.fields)
        locsout.(p.fields{k})=[];
    end
    indcombined=[];hroio=[];
    return
end

%obtain parameters
p=roiparser(varargin); 
if ischar(fields)
    p.fields={fields};
else
    p.fields=fields;
end

pd=p;

 %pass on selector
 if ~isempty(p.locselector)
     switch p.locselector.selector_filelist.selection
         case 'all'
             p.removeFilter=[p.removeFilter,{'filenumber'}];
         case 'layer'
         otherwise
             p.filenumber= p.locselector.selector_filelist.Value-2;
     end
     
     switch p.locselector.selector_filter.selection
         case 'layers'
             p.grouping='layer';
             p.layer=find(locData.getPar('sr_layerson'));
         case 'all u'
             p.grouping='ungrouped';
         case 'all g'
             p.grouping='grouped';
         otherwise %layers 1 to 3
             p.layer=p.locselector.selector_filter.Value-3;
     end     
     p.position=lower(p.locselector.selector_pos.selection);
 end


 %find implicit values, but ovrewrite with explicit definitions later
 fn=fieldnames(pd);
 for k=1:length(fn)
     if ~isempty(pd.(fn{k}))
         p.(fn{k})=pd.(fn{k});
     end 
 end
 
 %site as position
 if isa(p.position,'interfaces.SEsites')
     site=p.position;
    p.filenumber=site.info.filenumber;
    sr=locData.getPar('se_siteroi')/2;
    if isempty(sr)
        sr=300;
    end
    p.position=[site.pos(1:2) sr];
 end
 


%check if empty. necessary at beginning where xnm etc can be empty. Remove
%later?
if isempty(p.channel)&&strcmpi(p.grouping,'layer')&&isempty(p.layer)&&strcmpi(p.position,'default')&&isempty(p.filenumber)&&isempty(p.removeFilter)&&isempty(p.within)&&isempty(p.locselector)
locs=locData.loc;
    for k=1:length(p.fields)
        if isfield(locs,p.fields{k})
            locsout.(p.fields{k})=locs.(p.fields{k});
        elseif strcmp(p.fields{k},'ingrouped')
            locsout.(p.fields{k})=true(length(locData.grouploc.frame),1);
        elseif strcmp(p.fields{k},'inungrouped')
            locsout.(p.fields{k})=true(length(locData.loc.frame),1);
        else
            locsout.(p.fields{k})=[];
        end
    end
    hroio=[];
    indcombined=[];
    return
end

 hroio=[];
 
% remove filter in case field is defined
if ~iscell(p.removeFilter)
    p.removeFilter={p.removeFilter};
end
if isnumeric(p.position)||strcmp(p.position,'all')||strcmp(p.position,'default')
    p.removeFilter=[p.removeFilter,{'xnm','ynm'}];
end
if ~isempty(p.channel)
    p.removeFilter=[p.removeFilter,{'channel'}];
end
if ~isempty(p.filenumber)
    p.removeFilter=[p.removeFilter,{'filenumber'}];
end


%new approach:
% loop over layers
% for each layer get localizations, append them, irrespective of grouped,
% ungrouped
% define layer field
% set indgrouped, indungrouped: maybe matrix: N x layers
% if isempty(layer): all same, but no filtering. 
locstot=0;
layers=p.layer;
if isempty(layers)
    layers=1; %unfiltered locs
end

inungrouped=1;
ingrouped=1;


if any(strcmp(p.fields,'inlayerg')) || any(strcmp(p.fields,'ingrouped'))
    doinlayerg=true;
else
    doinlayerg=false;
end 
if any(strcmp(p.fields,'inlayeru')) || any(strcmp(p.fields,'inungrouped'))
    doinlayeru=true;
else
    doinlayeru=false;
end

for k=1:length(layers)
    switch p.grouping
        case {'layer',''}
            grouping=locData.isgrouped(layers(k));
        case 'grouped'
            grouping=true;
        case 'ungrouped'
            grouping=false;
 
    end
    if grouping && ~isempty(locData.grouploc)
        locs=locData.grouploc;
    else
        grouping=false;
        locs=locData.loc;
    end
    grouplayer(k)=grouping;

    if ~isempty(p.layer)%get filters
        if p.layer(k)>length(locData.layer)||p.layer(k)<1
            disp('getloc layer out of range')
            continue
        end
        filterold=locData.getFilter(p.layer(k),grouping); %use localizations if in any of the layers
        if ~isempty(p.removeFilter)
             filter=myrmfield(filterold,p.removeFilter);
             locData.setFilter(filter,p.layer(k),grouping);
        end
        indfilter=locData.inFilter(p.layer(k),grouping);
        locData.setFilter(filterold,p.layer(k),grouping);
    else
        indfilter=true(size(locs.xnm));
    end
    
     if ~isempty(p.channel)
         indfilterh=(locs.channel==p.channel(1));
         for c=2:length(p.channel)
            indfilterh=indfilterh|(locs.channel==p.channel(c));
         end
         indfilter=indfilter&indfilterh;
     end
     if ~isempty(p.filenumber)
         indfilterh=(locs.filenumber==p.filenumber(1));
         for f=2:length(p.filenumber)
            indfilterh=indfilterh|(locs.filenumber==p.filenumber(f));
         end
         indfilter=indfilter&indfilterh;
     end
     
    %get positions
     
    if isnumeric(p.position)
        pos=p.position(1:2);
        if length(p.position)==4
            sr_size=p.position(3:4)/2;
              indpos=abs(locs.xnm-pos(1))<sr_size(1) & abs(locs.ynm-pos(2))<sr_size(2);
              
        elseif length(p.position)==3
            indpos=(locs.xnm-pos(1)).^2+(locs.ynm-pos(2)).^2<=p.position(3)^2;

        else
            disp('wrong position parameter');
        end
    else
        switch p.position
            case {'all','default',''}
                indpos=true(size(locs.xnm));
            case 'roi'
                [indpos,hroio,strucout]=getinroi(locData,locs.xnm,locs.ynm,p.shiftxy);

                if isfield(strucout,'xnmline')
                    p.shiftxy(1:2)=0;
                   locs.xnmline=strucout.xnmline;
                   locs.ynmline=strucout.ynmline;
                end 
            case 'fov'
                pos=locData.getPar('sr_pos');
                sr_size=locData.getPar('sr_size');
                indpos=locs.xnm>pos(1)-sr_size(1) & locs.xnm<pos(1)+sr_size(1) & locs.ynm>pos(2)-sr_size(2) & locs.ynm<pos(2)+sr_size(2);
            otherwise %numerical position vector
                if isnumeric(p.position)
                    pos=p.position(1:2);
                    sr_size=p.position(3:4);
                    indpos=locs.xnm>pos(1)-sr_size(1) & locs.xnm<pos(1)+sr_size(1) & locs.ynm>pos(2)-sr_size(2) & locs.ynm<pos(2)+sr_size(2);
                else

                    disp('getlocs: position description not valid')
                end
        end

    end
    if ~isempty(p.within)
        indwithin=p.within;
        if length(indwithin)<length(indfilter)
            indwithin=grouped2ungrouped(locData,indwithin);
        elseif length(indwithin)>length(indfilter)
            indwithin=ungrouped2grouped(locData,indwithin);
        end
    else
        indwithin=true;
    end
    
    indcombined{k}=indfilter&indpos&indwithin;
    if doinlayerg
        inlayerg{k}=getindices(locData,indcombined{k},1);
        ingrouped=inlayerg{k}&ingrouped;
    end
    
    if doinlayeru
        inlayeru{k}=getindices(locData,indcombined{k},0);
        inungrouped=inlayeru{k}&inungrouped;
    end
    
    locslayer(k)=sum(indcombined{k});
    locstot=locstot+locslayer(k);

end 

if isempty(p.fields)||any(strcmpi(p.fields,'all'))
    p.fields=fieldnames(locs);
end

for k=1:length(p.fields)
    if isfield(locs,p.fields{k})
        locsout.(p.fields{k})=zeros(locstot,1,'like',locs.(p.fields{k}));
    else
         locsout.(p.fields{k})=zeros(locstot,1,'single');
    end
end

ind1=1;
for l=1:length(layers)
    if grouplayer(l)
        locs=locData.grouploc;
    else
        locs=locData.loc;
    end
    ind2=ind1+locslayer(l)-1;
    for k=1:length(p.fields)
         field=p.fields{k};
         if isfield(locs,field)
             vh=addshift(locs.(field),field,p.shiftxy);
             vh2=vh(indcombined{l});
            locsout.(field)(ind1:ind2)=vh2;         
         else
             locsout.(field)=[];
         end
    end 
    ind1=ind2+1;
end

if any(strcmp(p.fields,'ingrouped'))
    locsout.ingrouped=ingrouped;
end
if any(strcmp(p.fields,'inungrouped'))
    locsout.inungrouped=inungrouped;
end
if any(strcmp(p.fields,'indlayer'))
    locsout.indlayer=indcombined;
end
if any(strcmp(p.fields,'inlayerg'))
    locsout.inlayerg=inlayerg;
end
if any(strcmp(p.fields,'inlayeru'))
    locsout.inlayeru=inlayeru;
end
end

function out=addshift(in,field,shift)
out=in;
switch field
    case 'xnm'
        if shift(1)~=0
            out=in+shift(1);
        end
    case 'ynm'
        if shift(2)~=0
            out=in+shift(1);
        end
    case 'znm'
        if length(shift)>2&&shift(3)~=0
            out=in+shift(1);
        end      
end
end
function ind=getindices(obj,indcombined,isgrouped)
    if isgrouped
        fn=fieldnames(obj.grouploc);
        len=length(obj.grouploc.(fn{1}));
    else
        fn=fieldnames(obj.loc);
        len=length(obj.loc.(fn{1}));
    end
    if len==length(indcombined)
        ind=indcombined;
    elseif isgrouped
        ind=ungrouped2grouped(obj,indcombined);
    else
        ind=grouped2ungrouped(obj,indcombined);
    end   
end

function ind=ungrouped2grouped(obj,indcombined)
gind=unique(obj.loc.groupindex(indcombined));
ind=false(length(obj.grouploc.frame),1);
ind(gind)=true;
end

function ind=grouped2ungrouped(obj,indcombined)
ind=indcombined(obj.loc.groupindex);
end

function [indg,hroio,strucout]=getinroi(obj,vx,vy,shiftxy)
hroi=obj.getPar('sr_roihandle');
strucout=[];
if isa(hroi,'imroi')&&isvalid(hroi)  
        
    if isa(hroi,'imline')
        lw=obj.getPar('linewidth_roi')/2;
        pol=hroi.getPosition;
        dpol=pol(2,:)-pol(1,:);
        alpha=-atan2(dpol(1),dpol(2));

        len=sqrt(sum(dpol.^2))*1000/2;
        midp=mean(pol,1)*1000;
        [xr,yr]=rotcoord(vx-midp(1)+shiftxy(1),vy-midp(2)+shiftxy(2),alpha);
        indb=abs(xr)>lw|abs(yr)>len;
        indg=~indb;
        strucout.ynmline=xr;strucout.xnmline=yr;

    elseif isa(hroi,'impoint')
        lw=obj.parameters.linewidth_roi/2;
        pol=hroi.getPosition*1000;
        indg=abs(vx-pol(1))<lw&abs(vy-pol(2))<lw;
    else
        imbw=createMask(hroi,obj.getPar('sr_imagehandle'));
        sizeim=size(imbw);
        pos=obj.getPar('sr_pos');
        sizesr=obj.getPar('sr_size');
        pixrec=obj.getPar('sr_pixrec');
        xl=min(max(1,round((vx-pos(1)+sizesr(1))/pixrec)),sizeim(2));
        yl=min(max(1,round((vy-pos(2)+sizesr(2))/pixrec)),sizeim(1));
        linind=sub2ind(sizeim,yl,xl);
        indg=imbw(linind);
    end
    hroio=hroi;
else 
    srp=obj.getPar('sr_pos');
    ssize=obj.getPar('sr_size');
    pos=(srp(1:2)-ssize);
    
    hroio.getPosition=[pos(1) pos(2) 2*ssize(1) 2*ssize(2)]/1000;
    indg=vx>pos(1)&vx<pos(1)+2*ssize(1)&vy>pos(2)&vy<pos(2)+2*ssize(2);
    strucout=[];
end
end

function pres=roiparser(args)
p = inputParser;   
p.KeepUnmatched=true;
addParameter(p,'grouping','',@(x) any(myvalidatestring(x,{'grouped','ungrouped','layer'})));
addParameter(p,'layer',[],@isnumeric);
addParameter(p,'channel',[],@isnumeric);
addParameter(p,'filenumber',[],@isnumeric);
addParameter(p,'position','');
addParameter(p,'removeFilter',{});
addParameter(p,'within',[]);
addParameter(p,'shiftxy',[0,0,0]);
addParameter(p,'locselector',[]);
parse(p,args{:});
pres=p.Results;
if ~isempty(fieldnames(p.Unmatched))
    warning('locData.getlocs called with unidentified parameter - value pairs')
    disp(p.Unmatched)
end
end