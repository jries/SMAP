classdef Eval2loc<interfaces.DialogProcessor&interfaces.SEProcessor
%     Adds two field to the localization data containing the site number
%     and the cell number, respectively.
    methods
        function obj=Eval2loc(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
%             obj.inputParameters={'se_layerson'};
            obj.showresults=false;
        end
        
        function out=run(obj,p)  
            out=[];
            obj.setPar('undoModule','Eval2loc');
            notify(obj.P,'backup4undo');
            
            sites=obj.SE.sites;
            
            ld=obj.locData.loc;
            vals=0*ld.xnm+p.defaultval;
            sitenumbers=0*ld.xnm;
            roisize=obj.getPar('se_siteroi');
            fovsize=obj.getPar('se_sitefov');
            for k=1:length(sites)
                posxy=sites(k).pos(1:2);
                switch p.roimode.selection   
                    case 'site FoV'
                        position=[posxy fovsize fovsize];
                    case 'ROI square'
                        position=[posxy roisize roisize];
                    case 'ROI round'
                        position=[posxy roisize/2];
                end
                [l,ind]=obj.locData.getloc({'filenumber','xnm','ynm'},'Position',position);
                indh=ind & (obj.locData.loc.filenumber==sites(k).info.filenumber);
                assignedind=(sitenumbers>0) & indh;  %later populate manually if too slow
                
                site=sites(k);
                try
                    vhere=eval(p.property);
                catch err %did not work, not assigned
                    continue
                    disp([num2str(k) ' did not have property, ' p.property])
                end

                if any(assignedind)
                    newindall=find(indh);
                    oldsites=sitenumbers(assignedind);
                    usid=unique(oldsites);
                    posnh=sites(k).pos;
                    for kh=1:length(usid)
                        idh=obj.locData.SE.indexFromID(obj.locData.SE.sites,usid(kh));
                        posoldh=sites(idh).pos;
                       
                        dold=(l.xnm-posoldh(1)).^2+(l.ynm-posoldh(2)).^2;
                        dnew=(l.xnm-posnh(1)).^2+(l.ynm-posnh(2)).^2;
                        indold=dold<dnew;
                        newindall(indold)=0;
                    end
                    newindall(newindall==0)=[];
                    vals(newindall)=vhere;
                    sitenumbers(newindall)=sites(k).ID;
                else
                    vals(indh)=vhere;
                    sitenumbers(indh)=sites(k).ID;
                end
            end
            obj.locData.setloc(p.newfield,vals);
            obj.locData.regroup;
          
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end


function button_callback(a,b,obj)
site=obj.SE.currentsite;
ps.site.evaluation=site.evaluation;
ps.site.annotation=site.annotation;
ps.site.info=site.info;
ps.site.ID=site.ID;
ps.site.pos=site.pos;

str=browsefields(ps,'site');
 if ~isempty(str)
    val=eval(str);
    if isa(val,'double')&&length(val)>1
        numind=length(val(:));
        % for k=1:numind
        % numstr{k}=num2str(k);
        % end
        % answ=listdlg('ListString',numstr,'Name','Index','PromptString','select index','SelectionMode','single');
        % if ~isempty(answ)
        %     str=[str '(' num2str(answ) ')'];
        % end     
    end
    obj.guihandles.property.String=str;


    i1=strfind(str,'.');
    i2=strfind(str,'(');
    if isempty(i2)
        iend=length(str);
    else
        iend=i2(end)-1;
    end
    prop=str(i1(end)+1:iend);
    obj.guihandles.newfield.String=prop;

 end
end

function str=browsefields(prop,field)
if isempty(field)
    fn=fieldnames(prop);
elseif iscell(prop.(field))    
    v=prop.(field);
    if numel(v)==1
        index=1;
    else
        index=selectindexdialog(length(v));
    end
    fn=fieldnames(prop.(field){index});
    nextprop=prop.(field){index};
    field=[field '{' num2str(index) '}'];
elseif isstruct(prop.(field))
    fn=fieldnames(prop.(field));
    nextprop=prop.(field);
else
    str=field;
    return;
end
answ=listdlg('ListString',fn,'SelectionMode','single');
if isempty(answ)
    str='';
else
    field2=fn{answ};
    str=[field '.' browsefields(nextprop,field2)];
end

% if isstruc(prop)
% fn=fieldnames(prop)
% selectfield(fieldnames(obj.
end

function index=selectindexdialog(num)
for k=1:numind
   numstr{k}=num2str(k);
end
answ=listdlg('ListString',numstr,'Name','Index','PromptString','select index','SelectionMode','single');
if ~isempty(answ)
    str=[str '(' num2str(answ) ')'];
end     

end



function pard=guidef(obj)

pard.t1.object=struct('String','Write any evaluation or annotation to locData.loc','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=4;

pard.roimode.object=struct('String',{{'site FoV','ROI square','ROI round'}},'Style','popupmenu','Value',3);
pard.roimode.position=[2,1];
pard.roimode.Width=1;

pard.property.object=struct('String','','Style','edit');
pard.property.position=[3,1];
pard.property.Width=2;

pard.selectbutton.object=struct('String','select','Style','pushbutton','Callback',{{@button_callback,obj}});
pard.selectbutton.position=[3,3];

pard.newfieldt.object=struct('String','New field name:','Style','text');
pard.newfieldt.position=[4,1];
pard.newfieldt.Width=1;

pard.newfield.object=struct('String','','Style','edit');
pard.newfield.position=[4,2];
pard.newfield.Width=1;

pard.defaultvalt.object=struct('String','Not assigned:','Style','text');
pard.defaultvalt.position=[4,3];
pard.defaultvalt.Width=1;

pard.defaultval.object=struct('String','0','Style','edit');
pard.defaultval.position=[4,4];
pard.defaultval.Width=.5;


pard.plugininfo.description='Adds a field to the localization data containing the evaluation result.';
pard.plugininfo.type='ROI_Analyze';
end