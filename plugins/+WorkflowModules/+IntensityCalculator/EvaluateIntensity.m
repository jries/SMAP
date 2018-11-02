classdef EvaluateIntensity<interfaces.WorkflowModule
    %EvaluateIntensity Workflow module that calls intensity evaluators.
    properties
        evaluators
        intensities
        loccounter
        fields
        useevaluators
%         outputfilename
    end
    methods
        function obj=EvaluateIntensity(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
            obj.inputChannels=3;
        end
        function pard=guidef(obj)
             pard.plugininfo.type='WorkflowModule'; 
        end
        function makeGui(obj)
            makeGui@interfaces.WorkflowModule(obj);
            fw=obj.guiPar.FieldWidth;
            fh=obj.guiPar.FieldHeight;
            fs=obj.guiPar.fontsize;
            huitable=uitable('Units','pixels','Position',[fw,2*fh,fw,5*fh],'Parent',obj.handle,'FontSize',obj.guiPar.fontsize);
            huitable.ColumnEditable=[true false];
            huitable.RowName=[];
            huitable.ColumnName=[];
            huitable.ColumnWidth={fs,'auto'};
            huitable.CellSelectionCallback={@evalselect_callback,obj};
            huitable.Units='normalized';
            huitable.TooltipString='Select (check) evaluators to determine intensities';
            obj.guihandles.evalmodules=huitable;
            ev1={plugin('WorkflowModules','IntensityCalculator','roi2int_sumG')};
             obj.evaluators=ev1;
             ev2={plugin('WorkflowModules','IntensityCalculator','roi2int_fitG')};%roi2int_fitG;
             obj.evaluators(2)=ev2;
            p=obj.guiPar;
            p.Xpos=1;p.Vpos=1;p.Vrim=0;
            for k=1:length(obj.evaluators)
                ev=obj.evaluators{k};
                data{k,1}=true;
                data{k,2}=ev.info.name;
                hpanel=uipanel('Parent',obj.handle,'Units','pixels','Position',[2*fw,2*fh,2*fw,5*fh],'FontSize',fs,'Visible','off');
                 obj.children.(['panel_' num2str(k)])=ev;
                obj.guihandles.(['panel_' num2str(k)])=hpanel;
                ev.setGuiAppearence(p)
                ev.handle=hpanel;
                ev.makeGui;
                ev.handle.Units='normalized';
            end
            huitable.Data=data;
            obj.evaluators{1}.handle.Visible='on';
            obj.initGui;
        end
            
            
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.setInputChannels(3,'frame');  
        end
        function prerun(obj,data1,data2,data3)
            global EvaluateIntensity_intensity
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
                    for l=1:length(fields)
                    obj.fields={obj.fields{:} [fields{l} '2'] };
                    end
                    
                end
    
            end
               obj.fields={obj.fields{:} 'int_xpix', 'int_ypix', 'int_frame' };
            obj.intensities=single(0);
            EvaluateIntensity_intensity=single(0);
        end
        function dato=run(obj,data,p) %
            global EvaluateIntensity_intensity
            so1=2; %number of fields per evaluator
            so2=2;
            if ~isempty(data{1}.data)
                img=data{1}.data.img;
                bg=data{2}.data.img;
%                 img=d1.img;
%                 bg=d2.img;
                loc=data{3}.data;
                s=size(img);
                numl=s(3)/2;
                memincrease=1e5;  
                useevaluators=obj.useevaluators;
%                 intensities=obj.intensities;
                evaluators=obj.evaluators;
                loccounter=obj.loccounter;
                s=size(obj.intensities);
                if loccounter>s(1)-memincrease/4;
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
                        out2=evaluators{2}.evaluate(img,bg,loc.dx,loc.dy,loc.PSFxpix,loc.PSFypix);
                        EvaluateIntensity_intensity(loccounter+1:loccounter+numl,inds:inds+so2-1)=out2(1:numl,:);
                        inds=inds+so2;
                      EvaluateIntensity_intensity(loccounter+1:loccounter+numl,inds:inds+so2-1)=out2(numl+1:end,:);
                        inds=inds+so2;
                        
%                         intensities(loccounter+k,inds:inds+length(out)-1)=out;
%                         inds=inds+length(out);
%                         out=evaluators{2}.evaluate(im2,bg2,loc.dx(k+numl),loc.dy(k+numl),loc.PSFxpix(k+numl),loc.PSFypix(k+numl));
%                         intensities(loccounter+k,inds:inds+length(out)-1)=out;
                    end
                    EvaluateIntensity_intensity(loccounter+1:loccounter+numl,inds:inds+2)=horzcat(loc.x(1:numl),loc.y(1:numl),loc.frame(1:numl));
%                 end
                loccounter=loccounter+numl;
                
                obj.useevaluators=useevaluators;
%                 obj.intensities=intensities;
                obj.evaluators=evaluators;
                obj.loccounter=loccounter;
                dato=[];
                
            else
                dato=data{1};
%                 obj.output(data{1})
                if data{1}.eof
                    for k=1:length(obj.fields)
                        obj.locData.setloc(obj.fields{k},single(EvaluateIntensity_intensity(1:obj.loccounter,k)));
                    end
                    
                end
            end
        end

        function runold(obj,data)
%             disp(data{1}.frame)
            if ~isempty(data{1}.data)
                d1=data{1}.data;
                d2=data{2}.data;
                img=d1.img;
                bg=d2.img;
                loc=data{3}.data;
                s=size(img);
                numl=s(3)/2;
                memincrease=1e4;
                
                useevaluators=obj.useevaluators;
                intensities=obj.intensities;
                evaluators=obj.evaluators;
                loccounter=obj.loccounter;
                s=size(intensities);
                if loccounter>s(1)-memincrease/4;
                    intensities(loccounter+memincrease,length(obj.fields))=0;
                end
                
                for k=numl:-1:1
                    inds=1;
                    im1=img(:,:,k);
                    im2=img(:,:,k+numl);
                    bg1=bg(:,:,k);
                    bg2=bg(:,:,k+numl);
                    
                    if useevaluators(1)
                        out=evaluators{1}.evaluate(im1,bg1);
                        intensities(loccounter+k,inds:inds+length(out)-1)=out;
                        inds=inds+length(out);
                        out=evaluators{1}.evaluate(im2,bg2);
                        intensities(loccounter+k,inds:inds+length(out)-1)=out;
                        inds=inds+length(out);
                    end
                    if useevaluators(2)
                        out=evaluators{2}.evaluate(im1,bg1,loc.dx(k),loc.dy(k),loc.PSFxpix(k),loc.PSFypix(k));
                        intensities(loccounter+k,inds:inds+length(out)-1)=out;
                        inds=inds+length(out);
                        out=evaluators{2}.evaluate(im2,bg2,loc.dx(k+numl),loc.dy(k+numl),loc.PSFxpix(k+numl),loc.PSFypix(k+numl));
                        intensities(loccounter+k,inds:inds+length(out)-1)=out;
                    end
                end
                loccounter=loccounter+numl;
                
                obj.useevaluators=useevaluators;
                obj.intensities=intensities;
                obj.evaluators=evaluators;
                obj.loccounter=loccounter;
                
            else
                obj.output(data{1})
                if data{1}.eof
                    for k=1:length(obj.fields)
                        obj.locData.setloc(obj.fields{k},single(obj.intensities(1:obj.loccounter,k)));
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
