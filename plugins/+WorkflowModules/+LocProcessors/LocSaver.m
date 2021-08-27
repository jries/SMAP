classdef LocSaver<interfaces.WorkflowModule
%     Assembles a localization data structure from the fitted localizations
%     and saves it as a SMAP *.sml file. When fitting via a network,
%     fitting a local copy which is then moved to the destination can be
%     faster.
    properties
        timer
        filenumber
        fileinfo
        locDatatemp;
        deltaframes;
%         index;
%         numsaved
%         frames;
%         saveframes=100;
        savefields=struct('fieldnames',{{''}},'tosave',{{''}},'notsave',{{'PSFxerr','PSFyerr','bgerr','locpthompson','peakfindx','peakfindy'}});
        savefit
        
    end
    methods
       function obj=LocSaver(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
            obj.inputParameters={'loc_ROIsize'};
            obj.setInputChannels(1,[],'fitted localizations');
%             obj.propertiesToSave={'savefields'};
       end
       function savefit_callback(obj)
           savepar=obj.getPar('savefit');
           if isstruct(savepar)
               obj.savefit=copyfields(obj.savefit,savepar);
           elseif isnumeric(savepar) && savepar==0 %delete
               obj.savefit=[];
           end
           %delete command: 0: delete all.
           %otherwise: only structures are accepted and copied.
           
       end
        function pard=guidef(obj)
            pard.plugininfo.type='WorkflowModule'; 
            pard.plugininfo.description='Assembles a localization data structure from the fitted localizations and saves it as a SMAP *.sml file. When fitting via a network, fitting a local copy which is then moved to the destination can be faster.';
           
            pard.selectfields.object=struct('Style','pushbutton','String','Fields to save','Callback',{{@outputfields_callback,obj}});
            pard.selectfields.object.TooltipString='Select which fields to save. Use preview before.';
            pard.selectfields.position=[1,1];
            pard.selectfields.Width=1;
            pard.selectfields.Optional=true;
            
            pard.savelocal.object=struct('Style','checkbox','String','save local and copy','Value',0);
            pard.savelocal.object.TooltipString='Select this if you fit via a network and the saving of the localizations is very long (stauts bar stops for a long time at last frames).';
            pard.savelocal.position=[3,1];
            pard.savelocal.Width=2;
            pard.savelocal.Optional=true;
            
            pard.setoutputfile.object=struct('Style','pushbutton','String','Set:','Callback',{{@setoutputfile_callback,obj}});
            pard.setoutputfile.object.TooltipString='Set the output file';
            pard.setoutputfile.position=[2,1];
            pard.setoutputfile.Width=.5;
            pard.setoutputfile.Optional=true;
            
            pard.outputfile.object=struct('Style','edit','String','');
            pard.outputfile.object.TooltipString='output file name';
            pard.outputfile.position=[2,1.5];
            pard.outputfile.Width=1.5;
            pard.outputfile.Optional=true;
            
            pard.diffrawframest.object=struct('Style','text','String','Save every xxx raw frames');
            pard.diffrawframest.position=[4,1];
            pard.diffrawframest.Width=1.5;
            pard.diffrawframest.Optional=true;
            pard.diffrawframes.object=struct('Style','edit','String','500');
            pard.diffrawframes.position=[4,2.5];
            pard.diffrawframes.Width=0.5;
            pard.diffrawframes.Optional=true;
            
            pard.outputParameters={'diffrawframes'};
            pard.syncParameters={{'loc_outputfilename','outputfile',{'String','Value'}},{'savefit',[],{'String','Value'},@obj.savefit_callback}};
            
            
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)
            if obj.getPar('loc_preview')
                return
            end
%             numberofframes=obj.getPar('loc_fileinfo').numberOfFrames;
%             obj.deltaframes=floor(numberofframes/obj.saveframes);
%              obj.index=round(obj.deltaframes/2);
%             obj.numsaved=0;
%             obj.frames=struct('image',[],'frame',[]);
%             obj.frames=[];
            obj.locDatatemp=interfaces.LocalizationData;
            obj.locDatatemp.attachPar(obj.P);
            obj.locDatatemp.addfile;
            obj.filenumber=1;
             obj.locDatatemp.files.file=locsaveFileinfo(obj);  
             obj.locDatatemp.SE.files.info=obj.locDatatemp.files.file.info;
             obj.locDatatemp.SE.files.name=obj.locDatatemp.files.file.name;
             if contains(p.outputfile,'_sml.mat') %valid output name
                 obj.locDatatemp.files.file.name=p.outputfile;
             end
             p=obj.parent.getGuiParameters(true,true);
             p.name=obj.parent.pluginpath;
             
             obj.locDatatemp.addhistory(p);
             
             obj.fileinfo=obj.getPar('loc_cameraSettings');    %not used?  
              op=obj.getPar('overwrite_pixelsize');
              if ~isempty(op)
                  obj.fileinfo.cam_pixelsize_um=op;
              end
             
             obj.timer=tic;
             obj.run(); %clear variables;
        end
        function output=run(obj,data,p)
            persistent templocs numlocs
            if nargin<2
                templocs=[];numlocs=[];
                return
            end
            output=[];
            if obj.getPar('loc_preview')
                try 
                    obj.savefields.fieldnames=fieldnames(data.data); 
                end
                return
            end
            locs=data.data;%get;
            if ~isempty(locs)&&~isempty(locs.frame)
%                 maxfitdist=min(3.5,(p.loc_ROIsize-1)/2);
%                 indin=abs(locs.xpix-locs.peakfindx)<maxfitdist & abs(locs.ypix-locs.peakfindy)<maxfitdist;
                indin=true(size(locs.frame));
                if isfield(locs,'znm')
                    indin=indin&~isnan(locs.znm);
                end
                fn=fieldnames(locs);
                obj.savefields.fieldnames=fn;
                if isempty(templocs)
                    for k=1:length(fn)
                        templocs.(fn{k})(:,1)=locs.(fn{k})(indin);
                    end
                    numlocs=length(templocs.(fn{1}));
                else
                    newlocs=length(locs.(fn{1}));
                    if numlocs+newlocs>length(templocs.(fn{1}))
                        newlen=max(1000,2*(numlocs+length(locs.(fn{1}))));
                        for k=1:length(fn)
                            templocs.(fn{k})(newlen,1)=templocs.(fn{k})(end);
                        end
                    end
                    sindin=sum(indin);
%                     if sindin>0
                    for k=1:length(fn)
                        templocs.(fn{k})(numlocs+1:numlocs+sindin,1)=locs.(fn{k})(indin);
                    end
%                     end
                    numlocs=numlocs+sindin;
                end
%                 if locs.frame(end)>obj.index && obj.numsaved<obj.saveframes
%                     obj.numsaved=obj.numsaved+1;
%                     obj.index=obj.index+obj.deltaframes;
%                     if isempty(obj.frames)
%                         obj.frames=obj.getPar('loc_currentframe');
%                     else
%                         obj.frames(obj.numsaved)=obj.getPar('loc_currentframe');
%                     end
%                 end
            end
            
            
           
            
            if data.eof %save locs
                obj.status('saving localizations');drawnow
                if ~isempty(templocs)
                    locdat=interfaces.LocalizationData;
                    locdat.loc=fitloc2locdata(obj,templocs,1:numlocs);
                    obj.locDatatemp.addLocData(locdat);
                end
                filenameold=obj.locDatatemp.files.file(1).name;
                filename=filenameold;
                ind=2;
                while exist(filename,'file')
                    filename=[filenameold(1:end-7) num2str(ind) '_sml.mat'];
                    ind=ind+1;
                end
                mainfile=filename;
                if p.savelocal
                    filenameremote=filename;
                    filename=[pwd filesep 'temp.sml'];
                    mainfile=filename;
                end
                fitpar=obj.parent.getGuiParameters(true).children;
                fitpar.fittime=toc(obj.timer);
                fitpar.loadtifftime=obj.getPar('tiffloader_loadingtime');
                fitpar.processfittime=obj.getPar('tiffloader_fittime');
                fitpar.loc_globaltransform=obj.getPar('loc_globaltransform');
                obj.setPar('savefit',struct('fitparameters',fitpar)); obj.savefit_callback;
                try
                disp([num2str(length(obj.locDatatemp.loc.xnm)) ' localizations in ' num2str(fitpar.fittime) ' seconds.']);
                catch err
                    err
                end
                nosave=intersect(fieldnames(obj.locDatatemp.loc),obj.savefields.notsave);
                obj.locDatatemp.loc=rmfield(obj.locDatatemp.loc,nosave);
                
%                 average=obj.getPar('tiffloader_averagetiff');
                rawframes=obj.getPar('rawimagestack');
                if ~isempty(rawframes)
                    obj.locDatatemp.files.file.raw=rawframes;
                else
                    obj.locDatatemp.files.file.raw(1).image=0;
                    obj.locDatatemp.files.file.raw(1).frame=0;
                end
%                 obj.locDatatemp.files.file.raw(1).image=average;
%                 obj.locDatatemp.files.file.raw(1).frame=0;
                transformation=obj.getPar('loc_globaltransform');
                if isempty(transformation)
                    transformation=obj.getPar('loc_globaltransform3dcal'); %from spline fitter, only use if none selected in peak finder
                end
                if ~isempty(transformation)
                    obj.locDatatemp.files.file.transformation=transformation;
                end
                obj.locDatatemp.files.file.savefit=obj.savefit;
                
                if ~contains(filename,'nosave')
                try
                     obj.locDatatemp.savelocs(filename,[],struct('fitparameters',fitpar));
                catch err
                    [~,name,ext]=fileparts(filename);
                    filenamenew=[pwd filesep name ext];
                    obj.locDatatemp.savelocs(filenamenew,[],struct('fitparameters',fitpar));
                    warndlg('could not save sml file. Saved in local directory')
                    err
                end
                
                if p.savelocal
                    movefile(filename,filenameremote);
                end
                end
%               write to main GUI
%                 obj.locData.clear;
                obj.locData.setLocData(obj.locDatatemp);

                initGuiAfterLoad(obj);
                obj.setPar('mainfile',mainfile);
                [path,file]=fileparts(filename);
                try
                imageout=makeSRimge(obj,obj.locDatatemp);
                options.comp='jpeg';
                options.color=true;
                s=size(imageout);
                sr=ceil(s/16)*16;
                imageout(sr(1),sr(2),1)=0;
                saveastiff(uint16(imageout/max(imageout(:))*(2^16-1)),[path filesep file '.tif'],options)
                catch err
                    err
                end
                output=data;
            end
            
        end

    end
end

function outputfields_callback(a,b,obj)
tosave=obj.savefields.tosave;
notsave=obj.savefields.notsave;
fn=obj.savefields.fieldnames;
if isempty(fn)
    warndlg('please use preview before');
    return
end

fn{end+1}='filenumber';
on=num2cell(ones(length(fn),1));
loctest=cell2struct(on,fn,1);
fn=fieldnames(fitloc2locdata(struct('filenumber',1),loctest,1));

fn=setdiff(fn,{'xerrpix','yerrpix','PSFxpix','PSFypix','xpix','ypix','xnm','ynm','frame'});
tosave=setdiff(fn,notsave);
[save,inds]=checknames(fn,tosave);
if ~isempty(save)
    obj.savefields.tosave=tosave;
    obj.savefields.notsave=fn(~inds);
end
end

function imout=makeSRimge(obj,locDatatemp)
channelfile=[obj.getPar('SettingsDirectory')  '/workflows/FitTif_Channelsettings.mat'];
pall=load(channelfile);
p=pall.globalParameters;
p.lutinv=false;
p.sr_pixrec=20;
p.layer=1;
% p.gaussfac=0.4;
minx=min(locDatatemp.loc.xnm);maxx=max(locDatatemp.loc.xnm);
miny=min(locDatatemp.loc.ynm);maxy=max(locDatatemp.loc.ynm);
p.sr_pos=[(minx+maxx)/2 (miny+maxy)/2];
p.sr_size=[(maxx-minx)/2 (maxy-miny)/2];
p.sr_layerson=1;
rawimage=renderSMAP(locDatatemp.loc,p);
layers.images.finalImages=drawerSMAP(rawimage,p);
imoutt=displayerSMAP(layers,p);
imout=imoutt.image;
% imout=TotalRender(locDatatemp,pall.defaultChannelParameters);
end

function setoutputfile_callback(a,b,obj)
fn=obj.getSingleGuiParameter('outputfile');
if ~contains(fn,'sml.mat')
        fn='*sml.mat';
end
[f,p]=uiputfile(fn);
obj.setGuiParameters(struct('outputfile',[p f]));
    

end

