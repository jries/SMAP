classdef RoiAdder<interfaces.WorkflowModule
    properties
        maskrun
        mask
        preview
        
    end
    methods
        function obj=RoiAdder(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
            obj.propertiesToSave={'mask'};
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.guihandles.addroi.Callback={@addroi_callback,obj};
            obj.guihandles.removeroi.Callback={@removeroi_callback,obj};
            obj.guihandles.clearroi.Callback={@clearroi_callback,obj};
        end
        function prerun(obj,p)
            obj.preview=obj.getPar('loc_preview');
            if isempty(obj.mask)
%                 cf=obj.getPar('loc_currentframe');
%                 if isfield(cf,'image')
%                     sim=size(cf.image);
%                 else
                    cf=obj.getPar('loc_fileinfo');
                    if ~isempty(cf)&&myisfield(cf,'Width')
                        sim=[cf.Height cf.Width ];
                    else
                        sim=[0 0];
                    end
%                 end
                obj.maskrun=true(sim);  
            else
                obj.maskrun=obj.mask; 
            end
            if p.excluderim
            obj.setrim;
            end
 
        end
        function resetmask(obj,a,b)
%             cf=obj.getPar('loc_currentframe');
%             if isempty(cf)
                obj.mask=[];
%             else
%                 sim=size(cf.image);
%                 obj.mask=true(sim);   
%             end
        end
        
        function setrim(obj)
            dn=round((obj.getPar('loc_ROIsize')-1)/2);
            obj.maskrun(1:end,1:dn)=false;
            obj.maskrun(1:end,end-dn+1:end)=false;
            obj.maskrun(1:dn,1:end)=false;
            obj.maskrun(end-dn+1:end,1:end)=false;  
        end
        function dato=run(obj,data,p)

            img=data.data;%get;
            if isempty(img)
                dato=data;
                return;
            end
            if ~all(size(obj.maskrun)==size(img))
                obj.maskrun=imresize(obj.maskrun, size(img));
                %resize
                %set rim
                obj.setrim
            end
            
            img(~obj.maskrun)=-1;
           
            dato=data;%.copy;
            dato.data=img;%set(img);
            if obj.preview && data.frame==obj.getPar('loc_previewframe')
                figure(obj.getPar('loc_outputfig'));
                ax=gca;
                maxv=ax.CLim(2)/2;
                imagesc('CData',obj.maskrun*0+maxv,'AlphaData',double(~obj.maskrun)*0.5) 
            end  
        end
    end
end

function addroi_callback(a,b,obj)
p=obj.getGuiParameters;
obj.status('select roi and double clock on ROI when done...')
mask=getroi(obj.getPar('loc_outputfig'),p.roistyle.selection);
obj.status('ROI selected')
if isempty(obj.mask)||any(size(mask)~=size(obj.mask));
    obj.resetmask;
    obj.mask=mask;
else
    obj.mask=obj.mask|mask;
end
% figure(55);imagesc(obj.mask);
end


function removeroi_callback(a,b,obj)
p=obj.getAllParameters;
obj.status('select roi and double clock on ROI when done...')
mask=getroi(obj.getPar('loc_outputfig'),p.roistyle.selection);
obj.status('ROI selected')
if isempty(obj.mask)||any(size(mask)~=size(obj.mask));
    obj.resetmask;
    obj.mask=~mask;
else
    obj.mask=obj.mask&~mask;
end
end

function clearroi_callback(a,b,obj)
%     obj.mask=[];
    obj.resetmask;
end

function mask=getroi(figh,roistyle)
figure(figh);
ax=gca;
switch roistyle
    case 'rectangle'
        hroi=imrect(ax);
    case 'ellipse'
        hroi=imellipse(ax);
end
hpos=wait(hroi);
him=findobj(ax.Children,'Type','Image');
k=1;
while isempty(him(k).CData)
    k=k+1;
end

mask=createMask(hroi,him(k));
end

function pard=guidef
pard.addroi.object=struct('Style','pushbutton','String','ROI to include');
pard.addroi.position=[1,1];
pard.addroi.TooltipString=sprintf('Add a region of interest in which to fit. \n Use Preview before to display image in which to select ROI. \n After selcting ROI: double click on it to confirm.');
pard.removeroi.object=struct('Style','pushbutton','String','ROI to exclude');
pard.removeroi.position=[2,1];
pard.removeroi.TooltipString=sprintf('Add a region of interest in which NOT to fit. \n Use Preview before to display image in which to select ROI. \n After selcting ROI: double click on it to confirm.');

pard.clearroi.object=struct('Style','pushbutton','String','Clear ROIs');
pard.clearroi.position=[3,1];
pard.clearroi.TooltipString=sprintf('Clear all ROIs');
pard.roistyle.object=struct('Style','popupmenu','String','rectangle|ellipse');
pard.roistyle.position=[4,1];
pard.roistyle.TooltipString=sprintf('Use rectangular or elliptical ROI.');

pard.excluderim.object=struct('Style','checkbox','String','exclude rim','Value',1);
pard.excluderim.position=[5,1];
pard.excluderim.TooltipString=sprintf('Use only localizations sufficiently far rom edge of FoV.');
% pard.mask_store.object=struct('Style','text','String','');
pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Allows the user to select regions of interest on the preview image which to use for fitting, or which to exclude from fitting.';
end