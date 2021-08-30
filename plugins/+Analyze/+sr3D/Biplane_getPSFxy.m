classdef Biplane_getPSFxy<interfaces.DialogProcessor
%     associates size of the PSF in two channels to the localizations in
%     both channels
    methods
        function obj=Biplane_getPSFxy(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;   
           obj.showresults=false;
            obj.history=true;                
        end
        
        function out=run(obj,p)
           out=[];
            tload=load(p.Tfile);
           obj.setPar('undoModule','Biplane_get_PSFxy');
            notify(obj.P,'backup4undo');
            if ~isfield(tload,'transformation')
                out.error='selected transformation file does not have a valid transformation';
                return
            end
            
             p.datapart.selection='all';
             if p.allfiles
                 loc=obj.locData.loc;
                 keepold=true;
                 for k=1:length(obj.locData.files.file)
                     k
                     loc=get2ClocPSF(loc,tload.transformation,obj.locData.files.file(k),p,keepold);
                 end
             else
                 file=obj.locData.files.file(p.dataselect.Value);
                 keepold=false;
                loc=get2ClocPSF(obj.locData.loc,tload.transformation,file,p,keepold);
             end
            obj.locData.loc=copyfields(obj.locData.loc,loc);
            obj.locData.regroup;
             obj.setPar('locFields',fieldnames(obj.locData.loc))
        end
        
        function pard=guidef(obj)
            pard=guidef;
        end

        function initGui(obj)
            obj.addSynchronization('filelist_short',obj.guihandles.dataselect,'String')
            obj.addSynchronization('transformationfile',obj.guihandles.Tfile,'String');
        end

        function makeGui(obj)
            makeGui@interfaces.DialogProcessor(obj);
            obj.guihandles.loadbutton.Callback=@obj.loadbutton;
        end
        function loadbutton(obj,a,b)
            fn=obj.guihandles.Tfile.String;
            [f,path]=uigetfile(fn,'Select transformation file _T.mat');
            if f
                Tload=load([path f]);
                if ~isfield(Tload,'transformation')
                    msgbox('could not find transformation in file. Load other file?')
                end
                obj.guihandles.Tfile.String=[path f];
            end      
        end        
    end
end


function loco=get2ClocPSF(loc,transform,file,p,keepold)

p.datapart.selection='all (T->R)';
% p.datapart.selection='target';
infi=file.number==loc.filenumber;

inref=transform.getRef(loc.xnm,loc.ynm);
locref=copystructReduce(loc,infi&inref);
locth=copystructReduce(loc,infi&~inref);
loct=apply_transform_locs(locth,transform,file,p);

% only valid parts?
[iA,iB,uiA,uiB]=matchlocsall(renamefields(locref),renamefields(loct),0,0,1000);
length(iA)/(length(uiA)+length(uiB))
inffr=find(infi&inref);
infft=find(infi&~inref);


if ~keepold
    loco.PSFxnm_old=loc.PSFxnm;
    loco.PSFxnm=zeros(size(loc.xnm),'single');
    loco.PSFynm=zeros(size(loc.xnm),'single');
else
    loco=loc;
    if ~isfield(loco,'PSFxnm')
        loco.PSFxnm=zeros(size(loc.xnm),'single');
    end
    if ~isfield(loco,'PSFynm')
        loco.PSFynm=zeros(size(loc.xnm),'single');
    end
end


loco.PSFxnm(inffr(iA))=loc.PSFxnm(inffr(iA));
loco.PSFynm(inffr(iA))=loc.PSFxnm(infft(iB));
loco.PSFxnm(infft(iB))=loc.PSFxnm(inffr(iA));
loco.PSFynm(infft(iB))=loc.PSFxnm(infft(iB));
end

function loco=renamefields(loci)
loco.x=loci.xnm;
loco.y=loci.ynm;
loco.frame=loci.frame;

end

function pard=guidef
pard.texta.object=struct('String','dataset','Style','text');
pard.texta.position=[1,1];

pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[2,1];

pard.allfiles.object=struct('Style','checkbox','String','all files');
pard.allfiles.position=[2,2];
% pard.datapart.object=struct('Style','popupmenu','String','all|left|right|top|bottom');
% pard.datapart.position=[3,1];


pard.Tfile.object=struct('Style','edit','String','*_T.mat');
pard.Tfile.position=[8,1];
pard.Tfile.Width=3;

pard.loadbutton.object=struct('Style','pushbutton','String','load');
pard.loadbutton.position=[8,4];
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='associates size of the PSF in two channels to the localizations in both channels';
end