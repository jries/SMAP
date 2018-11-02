classdef Get2CIntImagesWF<interfaces.DialogProcessor
    % gets intensities from camera images at positions of localizations and
    % at transformed positions
    properties (Access=private)
        figure
    end
    methods
        function obj=Get2CIntImagesWF(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.history=true;
            obj.showresults=false;
        end
        
        function out=run(obj,p)
            out=[];
            if isempty(obj.figure)
                obj.figure=figure;
            end
            f=obj.figure;
            f.Visible='off';
            wffile='settings/workflows/get2CIntensityImagesWF.mat';
            wf=interfaces.Workflow(f,obj.P);
            wf.attachLocData(obj.locData);
            wf.makeGui;
            wf.load(wffile);
            
            transformation=loadtransformation(obj,p.Tfile,p.dataselect.Value);
            if isempty(transformation)
                out.error='selected transformation file does not have a valid transformation';
                return
            end
            obj.locData.files.file(p.dataselect.Value).transformation=transformation;
            file=obj.locData.files.file(p.dataselect.Value);
            fo=strrep(file.name,'_sml.mat','_dc_sml.mat');

            path=fileparts(file.name); %rather in top class, pass on
            %search for tiff images
            try
                if isfield(obj.locData.files.file(1).info,'imagefile')
                tiffile=obj.locData.files.file(1).info.imagefile;
                else
                    tiffile=obj.locData.files.file(1).info.basefile;
                end
                if ~exist(tiffile,'file')
                    disp('Get2CIntImagesWF ine 40: check if it works')
                    tiffile=strrep(tiffile,'\','/');
                    ind=strfind(tiffile,'/');
                    for k=1:length(ind)
                        tiffileh=[tiffile(1:ind(k)) '_b_' tiffile(ind(k)+1:end)];
                        if exist(tiffileh,'file')
                            tiffile=tiffileh;
                        end
                    end
                end
                if ~exist(tiffile,'file')
                    tiffile=[];
                end
            catch err
                tiffile=[];
            end
            
            if isempty(tiffile)
                [filen, path]=uigetfile([path filesep '*.tif'],file.name);
                tiffile=[path filen];
            end
%             tiffile='/Users/ries/Documents/Data/3Ddc/MTActin/02_MT680_phalloidin647_1/img_000039971_Default_000.tif';
            wf.module('TifLoader').addFile(tiffile);
            p.framestop=max(obj.locData.loc.frame)-1;
             wf.module('TifLoader').setGuiParameters(p);
             
            p.loc_blocksize_frames=p.filtert;
            p.loc_bg_dx=p.filterx;
            p.loc_subtractbg=true;
            wf.module('MedianBGcalculator').setGuiParameters(p);

            wf.module('IntLoc2pos').filestruc=file;
            wf.module('IntLoc2pos').setGuiParameters(p);

            pe=obj.children.evaluate.getGuiParameters(true);
            wf.module('EvaluateIntensity').setGuiParameters(pe,true);
            obj.setPar('loc_preview',false);
            wf.run;

            obj.locData.savelocs(fo);
            obj.locData.regroup;
            obj.setPar('locFields',fieldnames(obj.locData.loc))

        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            par.Vpos=3;
            par.Xpos=3;
            obj.children.evaluate=makeplugin(obj,{'WorkflowModules','IntensityCalculator','EvaluateIntensity'},obj.handle,par);
            obj.guihandles.loadbutton.Callback=@obj.loadbutton;

        end
        function loadbutton(obj,a,b)
            fn=obj.guihandles.Tfile.String;
            [f,path]=uigetfile(fn,'Select transformation file _T.mat');
            if f
                obj.guihandles.Tfile.String=[path f];
            end      
        end
    end
end




function pard=guidef

pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[1,1];
pard.dataselect.TooltipString=sprintf('Data set to be used');

pard.t1.object=struct('Style','text','String','Background filter');
pard.t1.position=[3,1];
pard.t1.Width=1;
pard.t1.TooltipString=sprintf('Median background filter. Select spatial and temporal step size');

pard.t2.object=struct('Style','text','String','dx');
pard.t2.position=[4,1];
pard.t2.Width=0.5;
pard.t2.TooltipString=sprintf('spatial step size (1-10, typical 3)');

pard.filterx.object=struct('Style','edit','String','3');
pard.filterx.position=[4,1.5];
pard.filterx.Width=0.5;
pard.filterx.TooltipString=pard.t2.TooltipString;

pard.t3.object=struct('Style','text','String','dt');
pard.t3.position=[5,1];
pard.t3.Width=0.5;
pard.t3.TooltipString=sprintf('temporal step size (e.g. 100)');

pard.filtert.object=struct('Style','edit','String','100');
pard.filtert.position=[5,1.5];
pard.filtert.Width=0.5;
pard.filtert.TooltipString=pard.t3.TooltipString;

pard.Tfile.object=struct('Style','edit','String','settings/temp/temp_T.mat');
pard.Tfile.position=[8,1];
pard.Tfile.Width=3;
pard.Tfile.TooltipString=sprintf('transformation file. Created e.g. with RegisterLocs');

pard.loadbutton.object=struct('Style','pushbutton','String','load');
pard.loadbutton.position=[8,4];
pard.loadbutton.TooltipString=pard.Tfile.TooltipString;

pard.syncParameters={{'filelist_short','dataselect',{'String'}},{'transformationfile','Tfile',{'String'}}};
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description=sprintf(['This plugin gets intensities from camera images at positions of localizations and at transformed positions \n',...
    'This plugin uses a transformation to find for every localization the position in the other channel and then determines the intensity in both channels.\n',...
    '1.	Load a transformation\n',...
    '2.	Per default, this plugin does median filtering. Select the spatial and temporal spacing for this (dx, dt).\n',...
    '3.	Select one or several plugins which determine the intensity:\n',...
    '\t a.	Roi2int_sum: uses a ROI (set size) to determine intensity, and a larger ROI for the background.\n',...
    '\t b.	Roi2int_fit: Uses a Gaussian fit to determine intensity and background. The position is fixed to the fitted position. You can use the fitted PSF size or fix it. If fit on BG is checked, the background is subtracted prior to fitting and the fit is performed with background set to zero. Otherwise the background is a fitting parameter.\n',...
    '4.	Press Run and when asked select the original raw camera images. The results are automatically saved with the _dc in the file name.\n']);
pard.plugininfo.name='2C intensities from images';
% pard.plugininfo.description
end



function plugino=makeplugin(obj,pluginpath,handle,pgui)
plugino=plugin(pluginpath{:});
plugino.attachPar(obj.P);
plugino.setGuiAppearence(pgui);
plugino.attachLocData(obj.locData);
plugino.handle=handle;
plugino.makeGui;
pluginname=pluginpath{end};
obj.children.(pluginname)=plugino;
end