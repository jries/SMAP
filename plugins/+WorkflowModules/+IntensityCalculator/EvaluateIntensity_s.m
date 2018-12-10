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
%         outputfilename
    end
    methods
        function obj=EvaluateIntensity_s(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
            obj.inputChannels=1;
        end
        function pard=guidef(obj)
             pard.plugininfo.type='WorkflowModule'; 
        end
        function makeGui(obj)
            makeGui@interfaces.WorkflowModule(obj);
            fw=obj.guiPar.FieldWidth;
            fh=obj.guiPar.FieldHeight;
            fs=obj.guiPar.fontsize;
            huitable=uitable('Units','pixels','Position',[5,2*fh,fw,5*fh],'Parent',obj.handle,'FontSize',obj.guiPar.fontsize);
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
%             ind=1;
            for k=1:length(obj.evaluators)
                if obj.makeevaluators(k)
                    ev=obj.evaluators{k};
                    data{k,1}=true;
                    data{k,2}=ev.info.name;
                    hpanel=uipanel('Parent',obj.handle,'Units','pixels','Position',[fw,2*fh,3*fw,5*fh],'FontSize',fs,'Visible','off');
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
            global EvaluateIntensity_intensity
            obj.extension=obj.getPar('intensity_channel');
            p=obj.getAllParameters;
            obj.useevaluators=[p.evalmodules.Data{:,1}];
            
            obj.loccounter=0;
            obj.fields={};
            for k=1:length(obj.evaluators)
                if obj.useevaluators(k)
                    obj.evaluators{k}.prerun;
                    fields=obj.evaluators{k}.info.fields;
                    for l=1:length(fields)
                    obj.fields={obj.fields{:} [fields{l} '1'] };
                    end
%                     for l=1:length(fields)
%                     obj.fields={obj.fields{:} [fields{l} '2'] };
%                     end
                    
                end
    
            end
               obj.fields={obj.fields{:} 'int_xpix', 'int_ypix', 'int_frame' ,'phot','bg'};
            obj.intensities=single(0);
            EvaluateIntensity_intensity=single(0);
        end
        function dato=run(obj,data,p) %
            global EvaluateIntensity_intensity
            so1=2; %number of fields per evaluator
            so2=2;
            so3=2;
            if ~isempty(data.data)
                img=data.data.img;
%                 bg=data{2}.data.img;
%                 img=d1.img;
%                 bg=d2.img;
                loc=data.data.info;
                s=size(img); 
                if length(s)==2
                    s(3)=1;
                end
                    
                numl=s(3);
                memincrease=1e4;  
                useevaluators=obj.useevaluators;
%                 intensities=obj.intensities;
                evaluators=obj.evaluators;
                loccounter=obj.loccounter;
                s=size(obj.intensities);
                s=size(EvaluateIntensity_intensity);
                if loccounter>s(1)-memincrease/4
                    EvaluateIntensity_intensity(loccounter+memincrease,length(obj.fields))=single(0);
                end
                
%                 for k=numl:-1:1
                    inds=1;
%                     im1=img(:,:,k);
%                     im2=img(:,:,k+numl);
%                     bg1=bg(:,:,k);
%                     bg2=bg(:,:,k+numl);
                    
                    if useevaluators(1)
                        out1=evaluators{1}.evaluate(img,bg);
                        out1p=out1(1:numl,:);
                        EvaluateIntensity_intensity(loccounter+1:loccounter+numl,inds:inds+so1-1)=out1p;
                         inds=inds+so1;
                         out1p=out1(numl+1:end,:);
                       EvaluateIntensity_intensity(loccounter+1:loccounter+numl,inds:inds+so1-1)=out1p;
                        inds=inds+so1;
                    end
                    if useevaluators(2)
                        out2=evaluators{2}.evaluate(img,loc.bg,loc.dx,loc.dy,loc.PSFxpix,loc.PSFypix);
                        EvaluateIntensity_intensity(loccounter+1:loccounter+numl,inds:inds+so2-1)=out2(1:numl,:);
                        inds=inds+so2;
%                       EvaluateIntensity_intensity(loccounter+1:loccounter+numl,inds:inds+so2-1)=out2(numl+1:end,:);
%                         inds=inds+so2;
                        
%                         intensities(loccounter+k,inds:inds+length(out)-1)=out;
%                         inds=inds+length(out);
%                         out=evaluators{2}.evaluate(im2,bg2,loc.dx(k+numl),loc.dy(k+numl),loc.PSFxpix(k+numl),loc.PSFypix(k+numl));
%                         intensities(loccounter+k,inds:inds+length(out)-1)=out;
                    end
                    if useevaluators(3)
                        out3=evaluators{3}.evaluate(img,loc);
                        EvaluateIntensity_intensity(loccounter+1:loccounter+numl,inds:inds+so3-1)=out3(1:numl,:);
                        inds=inds+so3;
                    end
                    EvaluateIntensity_intensity(loccounter+1:loccounter+numl,inds:inds+4)=horzcat(loc.x(1:numl),loc.y(1:numl),loc.frame(1:numl),loc.phot(1:numl),loc.bg(1:numl));
                    if isfield(loc,'groupindex')
                        EvaluateIntensity_intensity(loccounter+1:loccounter+numl,inds+5)=loc.groupindex(1:numl);
                    end
%                 end
                loccounter=loccounter+numl;
                
                
                obj.useevaluators=useevaluators;
%                 obj.intensities=intensities;
                obj.evaluators=evaluators;
                obj.loccounter=loccounter;
                 obj.setPar('fittedLocs',loccounter);
                dato=[];
                
            else
                dato=data;
%                 obj.output(data{1})
                if data.eof
                    groupindex=EvaluateIntensity_intensity(1:obj.loccounter,end);
                    [~,indsortg]=sort(groupindex);
                    Evis=single(EvaluateIntensity_intensity(indsortg,:));
                    for k=1:length(obj.fields)
                        % if grouped localizations: convert to ungrouped
                        % ones
                        if 0 %length(obj.locData.loc.xnm)>obj.loccounter
                            
                            %sort EvI by group index
                            %Evi(:,loc.gropuindex)
                            vh=Evis(obj.locData.loc.groupindex,k);
%                             ind=groupindex(obj.locData.loc.groupindex);
                        else
                            vh=single(EvaluateIntensity_intensity(1:obj.loccounter,k));
                        end
                        
                        obj.locData.setloc([obj.fields{k} obj.extension],vh);
                    end
                    
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
