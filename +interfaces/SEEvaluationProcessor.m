classdef SEEvaluationProcessor<interfaces.GuiModuleInterface & interfaces.LocDataInterface
    properties
        output
        name
        modulename
        number
        display
        site
         
    end

    methods
        function obj=SEEvaluationProcessor(varargin)
             obj@interfaces.GuiModuleInterface(varargin{:})
            obj.output.figure=-1;
            obj.outputParameters={'name','modulename'};
            
        end
        function redraw(obj)
            obj.site.image=[];
            obj.locData.SE.plotsite(obj.site)
        end
        
        function ax=setoutput(obj,name,clear)
            if nargin<3
                clear=false;
            end
            if obj.display
            if ~ishandle(obj.output.figure)
                obj.output.figure=figure;
                obj.output.tabgroup=uitabgroup('Parent',obj.output.figure);
                obj.output.tab=struct('none',0);
            end 
            
            tabnames=fieldnames(obj.output.tab);
            tabind=find(strcmpi(tabnames,name));
            if isempty(tabind)
                oldtab=obj.output.tabgroup.SelectedTab;
                tab=uitab(obj.output.tabgroup,'Title',name); 
                obj.output.tab.(name)=tab;
                delete(tab.Children);
                ax=axes('Parent',tab);
                if ~isempty(oldtab)
                obj.output.tabgroup.SelectedTab=oldtab;
                end
            end
            tab=obj.output.tab.(name);
            if clear
                delete(tab.Children);
                ax=axes('Parent',tab);
            end
            
            axs=findobj(tab.Children,'type','axes');
            ax=axs(1);
%             axes(ax);
            else
                h=figure(319);
                h.Visible='off';
                ax=axes('Parent',h);
            end
            if nargout==0
                axes(ax)
            end
        end
        
        function out=evaluate(obj,site)
            p=obj.getAllParameters;
            obj.site=site;
            out=obj.run(p);  
            site.evaluation.(obj.name)=out;    
            gp=rmfield(obj.getGuiParameters,{'classname','pluginpath'});
            site.evaluation.(obj.name).GuiParameters=gp; 
        end
        
        function [locsout,indloc]=getloc(obj,varargin) 
            [locsout,indloc]=obj.getLocs(varargin);
        end
        function [locsout,indloc]=getLocs(obj,varargin) 
            %uses locData.getloc
            %size: size of ROI. Otherwise it takes the whole FoV, not ROI
            %(from image rangex, rangey)
            %if one value: circular. vector: in x and y: rectangualr roi
            %size: freeroi: defined roi, or if not defined then site_roi
            %(circular)
            p=roiparser(varargin); 
            fields=varargin{1};
            fields={fields{:} ,'xnm','ynm','filenumber'};
            parameters=varargin(2:end);
            inds=find(strcmp(parameters,'size'));
            if ~isempty(inds)
                parameters(inds:inds+1)=[];
            end
            
                sx=obj.site.image.rangex;
                fovsize=(sx(2)-sx(1))*1000*[1,1];             
            if isempty(p.position) %no position specified, that is the usual case

                pos=obj.site.pos;
                if ~isempty(p.size)&&isnumeric(p.size)
                    if length(p.size)==1
                        fovsize=p.size(1);
                    else
                    fovsize=p.size;
                    end
                end
                posgetloc=[pos(1), pos(2), fovsize];
                parameters=[parameters {'position', posgetloc}];
            else
                pos=p.position;
            end
            [locsh, indloc]=obj.locData.getloc(fields,parameters{:},'removeFilter',{'filenumber'});
        
            inroi=true;
            if ischar(p.size)&&contains(p.size,'freeroi')
                hroi=obj.locData.SE.processors.preview.hlines.line3;
                if ~isempty(hroi)&&isvalid(hroi)
                    inroi=hroi.inROI(double(locsh.xnm/1000),double(locsh.ynm/1000));
                    p.size=inf;
                else
                    p.size=obj.getPar('se_siteroi')/2;  
                end
            end
            
            locsh.xnmr=locsh.xnm-pos(1);
            locsh.ynmr=locsh.ynm-pos(2);
            
            
            
            [locsh.xnmrot,locsh.ynmrot]=rotcoord(locsh.xnmr,locsh.ynmr,obj.site.annotation.rotationpos.angle/180*pi);
            
            if length(p.size)==1 %radius defined
                indroi=locsh.xnmr.^2+locsh.ynmr.^2<=p.size.^2;
            else
                if isempty(p.size)
                    p.size=[fovsize,fovsize];
                end
                indroi=locsh.xnmrot>-p.size(1) & locsh.xnmrot<p.size(1)&locsh.ynmrot>-p.size(2) & locsh.ynmrot<p.size(2);
            end
            indfile=locsh.filenumber==obj.site.info.filenumber;
            indgood=indfile&indroi&inroi;
            if iscell(indloc)
                indout=[];
                disp('not implemented: SEEvaluation procesor line 142')
            else
            indloc(indloc)=indgood;
            end
            
            fn=setdiff(fieldnames(locsh),{'ingrouped','inungrouped'});
            
             for k=1:length(fn)
                 field=fn{k};
                 if ~isempty(locsh.(field))
%                  if isfield(locs,field)
                    locsout.(field)=locsh.(field)(indgood);
                 else
                      locsout.(field)=[];
                 end
%                  end
             end  
             if isfield(locsh,'ingrouped')
                 locsh.ingrouped(locsh.ingrouped)=indgood;
             end
             if isfield(locsh,'inungrouped')
                 locsh.inungrouped(locsh.inungrouped)=indgood;
             end             
             locsout=copyfields(locsout,locsh,{'ingrouped','inungrouped'});
        end
    end
end

function pres=roiparser(args)

p = inputParser; 
p.KeepUnmatched=true;
addOptional(p,'fields',{});
%  addOptional(p,'fields',{},@(x) any(validatestring(x,{'cell'})));
addParameter(p,'grouping','ungrouped',@(x) any(validatestring(x,{'grouped','ungrouped'})));
addParameter(p,'layer',[]);
addParameter(p,'channel',[]);
addParameter(p,'rotate',false);
addParameter(p,'size',[]);
addParameter(p,'position',[]);
   parse(p,args{:});
   pres=p.Results;
end

