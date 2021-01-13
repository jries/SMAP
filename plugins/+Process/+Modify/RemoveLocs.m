classdef RemoveLocs<interfaces.DialogProcessor
%     remove localizations inside / outside a user-defined ROI
    methods
        function obj=RemoveLocs(varargin)   
            obj@interfaces.DialogProcessor(varargin{:});  
        end
        
        function out=run(obj,p)
            obj.setPar('undoModule','RemoveLocs');
            notify(obj.P,'backup4undo');
            [locs,indroi]=obj.locData.getloc('filenumber','position','roi','grouping','ungrouped');
            layers=find(obj.getPar('sr_layerson'));
%             finroi=find(indroi);
            if p.allfiles
%                 infile=true(size(locs.filenumber));
                infileall=true(size(obj.locData.loc.filenumber));
            else
%                 infile=false(size(locs.filenumber));
                infileall=false(size(obj.locData.loc.filenumber));
                for k=1:length(layers)
                    pl=obj.getPar(['layer' num2str(layers(k)) '_']);
%                     infile=infile|locs.filenumber==pl.ch_filelist.Value;
                    infileall=infileall|obj.locData.loc.filenumber==pl.ch_filelist.Value;
                end
            end
            
            
            switch p.roic.Value
            
                case 1
                    indout=indroi;
                    indout((~infileall))=false;
%                     indout=(indroi & infile);
                case 2
                    indout=~indroi;
                    indout((~infileall))=false;
%                     indout=(~indroi & infile);
                case 3
                    if p.allfiles
                        [~,indroi]=obj.locData.getloc('xnm','position','roi','layer',layers,'grouping','ungrouped','removeFilter',{'filenumber'});
                    else
                        [~,indroi]=obj.locData.getloc('xnm','position','roi','layer',layers,'grouping','ungrouped');
                    end
                    indout=~indroi;      
            end
            if p.setprop
                if isfield(obj.locData.loc,'active')
                    obj.locData.setloc('active',obj.locData.loc.active & ~indout);
                else
                   obj.locData.setloc('active', ~indout); 
                end
            else
                obj.locData.removelocs(indout);
            end
            
            obj.locData.regroup;   
            out=[];
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function updateGui(obj,event,data)
            if ~isempty(obj.locData.loc)
            fn=fieldnames(obj.locData.loc);
            obj.guihandles.fieldselect.String=fn;
            end
        end       

    end
end




function pard=guidef
pard.roic.object=struct('String',{{'remove inside ROI','remove outside ROI','keep visible inside ROI'}},'Style','popupmenu');
pard.roic.position=[1,1];
pard.roic.Width=2;
pard.roic.object.TooltipString='';

pard.allfiles.object=struct('String','all files','Style','checkbox');
pard.allfiles.position=[2,1];
pard.allfiles.Width=2;

pard.setprop.object=struct('String','set property: active','Style','checkbox');
pard.setprop.position=[3,1];
pard.setprop.Width=2;


pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='remove localizations inside / outside a user-defined ROI';
end