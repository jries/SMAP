classdef EvaluateIntensity_s<interfaces.WorkflowModule
    %EvaluateIntensity Workflow module that calls intensity evaluators.
    properties
        evaluators
        intensities
        loccounter
        fields
        useevaluators
        extension
        makeevaluators=true(3,1);
        peval
        EvaluateIntensity_intensity
        timershow
    end
    methods
        function obj=EvaluateIntensity_s(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
            obj.inputChannels=1;
        end
        function pard=guidef(obj)
             pard.plugininfo.type='WorkflowModule'; 
             pard.plugininfo.description='EvaluateIntensity Workflow module that calls intensity evaluators.';
        end
        function makeGui(obj)
            makeGui@interfaces.WorkflowModule(obj);
            fw=obj.guiPar.FieldWidth;
            fh=obj.guiPar.FieldHeight;
            fs=obj.guiPar.fontsize;
            huitable=uitable('Units','pixels','Position',[5,0,fw,4.5*fh],'Parent',obj.handle,'FontSize',obj.guiPar.fontsize);
            huitable.ColumnEditable=[true false];
            huitable.RowName=[];
            huitable.ColumnName=[];
            huitable.ColumnWidth={fs,fw-fs-5};
            huitable.CellSelectionCallback={@evalselect_callback,obj};
            huitable.Units='normalized';
            huitable.TooltipString='Select (check) evaluators to determine intensities';
            obj.guihandles.evalmodules=huitable;
            ev1={plugin('WorkflowModules','IntensityCalculator','roi2int_sumG')};
             obj.evaluators=ev1;
             ev2={plugin('WorkflowModules','IntensityCalculator','roi2int_fitG')};%roi2int_fitG;
             obj.evaluators(2)=ev2;
            ev3={plugin('WorkflowModules','IntensityCalculator','roi2int_expPSF')};%roi2int_fitG;
             obj.evaluators(3)=ev3;
            p=obj.guiPar;
            p.Xpos=1;p.Vpos=1;p.Vrim=0;
            for k=1:length(obj.evaluators)
                if obj.makeevaluators(k)
                    ev=obj.evaluators{k};
                    data{k,1}=true;
                    data{k,2}=ev.info.name;
                    hpanel=uipanel('Parent',obj.handle,'Units','pixels','Position',[fw,0,3*fw,4.5*fh],'FontSize',fs,'Visible','off');
                    obj.children.(['panel_' num2str(k)])=ev;
                    obj.guihandles.(['panel_' num2str(k)])=hpanel;
                    ev.setGuiAppearence(p)
                    ev.handle=hpanel;
                    ev.attachPar(obj.P);
                    ev.makeGui;
                    ev.handle.Units='normalized';
                end
            end
            huitable.Data=data;
            obj.evaluators{find(obj.makeevaluators,1,'first')}.handle.Visible='on';
            obj.initGui;
        end
            
            
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.setInputChannels(1,'frame');  
        end
        function prerun(obj,data1,data2,data3)
%             global EvaluateIntensity_intensity timershow
            extension=obj.getPar('intensity_channel');
            if ~isempty(extension)
                obj.extension=extension;
            end
            p=obj.getAllParameters;
            for k=1:size(p.evalmodules.Data,1)
                if isempty(p.evalmodules.Data{k,1})
                    obj.useevaluators(k)=0;
                else
                    obj.useevaluators(k)=p.evalmodules.Data{k,1};
                end
            end
            
            obj.loccounter=0;
            obj.fields={};
            for k=1:length(obj.evaluators)
                obj.peval{k}=obj.evaluators{k}.getAllParameters;
                if obj.useevaluators(k)
                    obj.evaluators{k}.extension=obj.extension;
                    obj.evaluators{k}.prerun(obj.peval{k});
                    fields=obj.evaluators{k}.info.fields;
                    for l=1:length(fields)
                        obj.fields={obj.fields{:} [fields{l}] };
                    end 
                end
                
            end
                obj.fields={obj.fields{:} 'int_xpix' 'index' };
            obj.intensities=single(0);
            obj.EvaluateIntensity_intensity=single(0);
            obj.timershow=tic;
        end
        function dato=run(obj,data,p) %
%             global EvaluateIntensity_intensity timershow
            so=2;
%             EvaluateIntensity_intensity=obj.EvaluateIntensity_intensity;
            if ~isempty(data.data)
                img=data.data.img;
                loc=data.data.info;
                s=size(img); 
                if length(s)==2
                    s(3)=1;
                end
                    
                numl=s(3);
                memincrease=1e7;  
                useevaluators=obj.useevaluators;
                evaluators=obj.evaluators;
                loccounter=obj.loccounter;
                
                s=size(obj.EvaluateIntensity_intensity);
                if loccounter>s(1)-memincrease/4
                    obj.EvaluateIntensity_intensity(loccounter+memincrease,length(obj.fields))=single(0);
                end
                
                inds=1; 
                for ev=1:length(evaluators)
                    if useevaluators(ev)
                        out=evaluators{ev}.evaluate(obj.peval{ev},img,loc);
                       obj.EvaluateIntensity_intensity(loccounter+1:loccounter+numl,inds:inds+so-1)=out(1:numl,:);
                        inds=inds+so;
                    end
                end
                obj.EvaluateIntensity_intensity(loccounter+1:loccounter+numl,inds+0)=single(loc.xpix(1:numl)); 
                obj.EvaluateIntensity_intensity(loccounter+1:loccounter+numl,inds+1)=single(loc.ind(1:numl));                     
                if isfield(loc,'groupindex')
                    obj.EvaluateIntensity_intensity(loccounter+1:loccounter+numl,inds+2)=loc.groupindex(1:numl);
                end
                loccounter=loccounter+numl;
                           
                obj.useevaluators=useevaluators;
                obj.evaluators=evaluators;
                obj.loccounter=loccounter;
                if toc(obj.timershow)>1
                    obj.setPar('fittedLocs',loccounter);
                    obj.timershow=tic;
                end
                dato=[];
%                 obj.EvaluateIntensity_intensity=EvaluateIntensity_intensity;
            else
                dato=data;
            end
                if data.eof
                    groupindex=obj.EvaluateIntensity_intensity(1:obj.loccounter,end);
                    [~,indsortg]=sort(groupindex);
                    Evis=single(obj.EvaluateIntensity_intensity(indsortg,:));
                    for k=1:length(obj.fields)
                        % if grouped localizations: convert to ungrouped
                        % ones
                        if 0 %length(obj.locData.loc.xnm)>obj.loccounter
                            
                            %sort EvI by group index
                            %Evi(:,loc.gropuindex)
                            vh=Evis(obj.locData.loc.groupindex,k);
%                             ind=groupindex(obj.locData.loc.groupindex);
                        else
                            vh=single(obj.EvaluateIntensity_intensity(1:obj.loccounter,k));
                        end
                        
                        obj.locData.setloc([obj.fields{k} obj.extension],vh);
                    end
                    
                end
            
        end

    end
end

function evalselect_callback(a,data,obj)
row=data.Indices(1);
for k=1:length(obj.evaluators)
    obj.evaluators{k}.handle.Visible='off';
end
obj.evaluators{row}.handle.Visible='on';

end
