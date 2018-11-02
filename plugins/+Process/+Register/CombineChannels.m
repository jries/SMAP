classdef CombineChannels<interfaces.DialogProcessor
    % CombineChannels applies transformation to localizations and combines
    % associated localizations into one
    methods
        function obj=CombineChannels(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
%             obj.inputParameters={'cam_pixelsize_nm'};
            obj.showresults=false;
            obj.history=true;
        end
        
        function out=run(obj,p) 
            out=[];
            obj.setPar('undoModule','CombineChannels');
            notify(obj.P,'backup4undo');
            load(p.Tfile)
            if ~exist('transformation','var')
                out.error='selected transformation file does not have a valid transformation';
                return
            end
                
            if p.allfiles
                filenumbers=1:length(obj.locData.files.file);
            else
                filenumbers=p.dataselect.Value;
            end
            obj.locData.addloc('matched');
            for f=1:length(filenumbers)
                file=obj.locData.files.file(filenumbers(f));
                locDatar=obj.locData.copy;
                thisfile=locDatar.loc.filenumber==filenumbers(f);
                locDatar.removelocs(~thisfile);
%                 locDatat=locDatar.copy;
                
                indref=transformation.getRef(locDatar.loc.xnm,locDatar.loc.ynm);
                loc=copystructReduce(locDatar.loc,indref);
                loct=copystructReduce(locDatar.loc,~indref);
                loctt=apply_transform_locs(loct,transformation,file);%,struct('datapart',struct('selection','target')));
                [iA,iB,uiA,uiB]=matchlocsall(renamefields(loc),renamefields(loctt),0,0,500);
                loctt=copyfields(loct,loctt,fieldnames(loc));
                matchedA=copystructReduce(loc,iA);
                matchedB=copystructReduce(loctt,iB);
                unmatchedA=copystructReduce(loc,uiA);
                unmatchedB=copystructReduce(loctt,uiB);
                fn=fieldnames(loc);
                if isfield(matchedA,'locprecnm')
                    wA=1./matchedA.locprecnm;
                    wB=1./matchedB.locprecnm;
                else
                    photA=matchedA.phot;
                    photB=matchedB.phot;
                    wA=sqrt(photA);
                    wB=sqrt(photB);
                end
                for k=1:length(fn)
                    switch fn{k}
                        case {'xnm','ynm','znm','PSFxnm','PSFynm','locprecnm','locprecznm'}
                            av=(matchedA.(fn{k}).*wA+matchedB.(fn{k}).*wB)./(wA+wB);
                        case {'phot','bg'}
                            av=(matchedA.(fn{k})+matchedB.(fn{k}));
                        otherwise
                            av=matchedA.(fn{k});
                    end

                    locout.(fn{k})=vertcat(av,unmatchedA.(fn{k}),unmatchedB.(fn{k}));
                end
                locout.matched=vertcat(ones(size(matchedA.(fn{1}))),zeros(size(unmatchedA.(fn{1}))),zeros(size(unmatchedB.(fn{1}))));

%                 obj.locData.loc=copyfields(obj.locData.loc,loc);
            
%             fn=fieldnames(obj.locData.loc);
%             for k=1:length(fn)
%                 obj.locData.loc.(fn{k})(thisfile)=locout.(fn{k});
%             end
%             end
            obj.locData.removelocs(thisfile);
            obj.locData.addLocData(locout);
            end
            obj.locData.regroup;
            obj.setPar('locFields',fieldnames(obj.locData.loc));
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function makeGui(obj)
            makeGui@interfaces.DialogProcessor(obj);
        end
        function initGui(obj)
            obj.addSynchronization('transformationfile',obj.guihandles.Tfile,'String');
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
                obj.setPar('transformationfile',[path f]);
            end      
        end  
    end
end

function loco=renamefields(loci)
loco.x=loci.xnm;
loco.y=loci.ynm;
loco.frame=loci.frame;
end


function outind=inref(x,y,separator)
inind=x>0&x<separator(1)&y>0&y<separator(2);
outind=~inind;
end


function pard=guidef(obj)
pard.texta.object=struct('String','dataset:','Style','text');
pard.texta.position=[1,1];

pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[1,1.7];
pard.dataselect.Width=2;

pard.allfiles.object=struct('Style','checkbox','String','transform all files','Value',0);
pard.allfiles.position=[1,3.7];
pard.allfiles.Width=1.3;
% 
% pard.textb.object=struct('String','transform:','Style','text');
% pard.textb.position=[2,1];

% pard.datapart.object=struct('Style','popupmenu','String',{{'all (T->R)','all (R->T)','reference','target'}});
% pard.datapart.position=[2,1.7];
% 


pard.Tfile.object=struct('Style','edit','String','settings/temp/temp_T.mat');
pard.Tfile.position=[7,1];
pard.Tfile.Width=3;

pard.loadbutton.object=struct('Style','pushbutton','String','load T','Callback',@obj.loadbutton);
pard.loadbutton.position=[7,4];

pard.syncParameters={{'filelist_short','dataselect',{'String'}}};

pard.plugininfo.description='ApplyTransform applies transformation to localizations or associated Tif images';
pard.plugininfo.name='Combine Channels';
pard.plugininfo.type='ProcessorPlugin';
end