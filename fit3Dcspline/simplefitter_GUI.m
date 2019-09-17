%  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
%  Heidelberg.
%  
%  This file is part of GPUmleFit_LM Fitter.
%  
%  GPUmleFit_LM Fitter is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  
%  GPUmleFit_LM Fitter is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with GPUmleFit_LM Fitter.  If not, see <http://www.gnu.org/licenses/>.
%  
%  
%  Additional permission under GNU GPL version 3 section 7

%%
classdef simplefitter_GUI<handle
    properties
        guihandles
        smappos
        mij
    end
    methods
        function obj=simplefitter_GUI(varargin)  
            %constructur: make GUI   
            if ~isdeployed
            addpath('shared')
            addpath('bfmatlab')
            end
            javaaddpath('ImageJ/plugins/bioformats_package.jar')
            
            figureheight=630;
            h=figure('Name','3D fitter','MenuBar','none','ToolBar','none');
            initPosition = h.Position;
            h.Position=[initPosition(1), initPosition(2)- figureheight+initPosition(4),450, figureheight];
            top=h.Position(4)-10;
            vsep=24;
            if ispc
                fontsize=12;
            else 
                fontsize=14;
            end
            xpos1=10;
            xw=105;
            hatitle='left';
            obj.guihandles.handle=h;
            obj.guihandles.title=uicontrol('style','text','String','Fit image stack with experimental PSF model. (c) 2017 Ries lab','Position',[xpos1,top-vsep+10,xw*4.5,vsep],'FontSize',10,'HorizontalAlignment',hatitle,'FontWeight','bold');
            
            if ismac
                obj.guihandles.mac=uicontrol('style','text','String','GPU fit only works on Windows. On Macs the CPU fitter is used.','Position',[xpos1,top-vsep-5,xw*4.5,vsep],'FontSize',10,'HorizontalAlignment',hatitle,'FontWeight','bold');
            end
   
             obj.guihandles.filetxt=uicontrol('style','text','String','Data and calibration files:','Position',[xpos1,top-2*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment',hatitle,'FontWeight','bold');
            obj.guihandles.loadert=uicontrol('style','text','String','Loader','Position',[xpos1,top-3*vsep,xw*1.5,vsep],'FontSize',fontsize);
            obj.guihandles.loader=uicontrol('style','popupmenu','String',{'simple tif','ome loader','ImageJ'},'Position',[xpos1+1.5*xw,top-3*vsep,xw*2.5,vsep],'FontSize',fontsize,'Callback',@obj.changeloader_callback,'Value',3);
                                  
            obj.guihandles.selectfiles=uicontrol('style','pushbutton','String','Load raw images','Position',[xpos1,top-4*vsep,xw*2,vsep],'FontSize',fontsize,'Callback',{@obj.selectfiles_callback,1});
             obj.guihandles.selectfiles.TooltipString='Select raw image file (tiff stack).';
            obj.guihandles.imagefile=uicontrol('style','edit','String','','Position',[xpos1+2*xw,top-4*vsep,xw*2,vsep],'FontSize',fontsize);
            obj.guihandles.imagefile.TooltipString='Tiff stack containing the single-molecule data.';
            
            obj.guihandles.selectcoeff=uicontrol('style','pushbutton','String','Load calibration','Position',[xpos1,top-5*vsep,xw*2,vsep],'FontSize',fontsize,'Callback',{@obj.selectfiles_callback,2});
            obj.guihandles.selectcoeff.TooltipString='Select 3D calibration file generated with the tool calibrat3D_GUI.';
            obj.guihandles.calfile=uicontrol('style','edit','String','','Position',[xpos1+2*xw,top-5*vsep,xw*2,vsep],'FontSize',fontsize);
            obj.guihandles.calfile.TooltipString='File containing the 3D calibration. *_3dcal.mat. Leave empty for Gaussian fitting.';
             
            obj.guihandles.isscmos=uicontrol('style','checkbox','String','sCMOS','Position',[xpos1,top-6*vsep,xw,vsep],'FontSize',fontsize);   
            obj.guihandles.selectscmos=uicontrol('style','pushbutton','String','Load var map','Position',[xpos1+xw,top-6*vsep,xw*1,vsep],'FontSize',fontsize,'Callback',{@obj.selectfiles_callback,3});   
            obj.guihandles.selectscmos.TooltipString='Select sCMOS variance map (in ADU^2) of same size ROI on chip as image stack';
            obj.guihandles.scmosfile=uicontrol('style','edit','String','','Position',[xpos1+2*xw,top-6*vsep,xw*2,vsep],'FontSize',fontsize);
            obj.guihandles.scmosfile.TooltipString='Tiff/.mat image containing sCMOS variance map (same ROI on camera as tiff).';
            


            obj.guihandles.camtext=uicontrol('style','text','String','Acquisition parameters:','Position',[xpos1,top-8*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment',hatitle,'FontWeight','bold');
            
            ha='right';
            
            obj.guihandles.conversiont=uicontrol('style','text','String','conversion (e-/ADU)','Position',[xpos1,top-9*vsep,xw*1.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.conversion=uicontrol('style','edit','String','0.1','Position',[xpos1+xw*1.5,top-9*vsep,xw*.5,vsep],'FontSize',fontsize);
            obj.guihandles.conversion.TooltipString=sprintf('conversion factor = conv/EMgain. conv is the gain stated in the camera spec sheet (e-/ADU)');

            obj.guihandles.conversiont.TooltipString=obj.guihandles.conversion.TooltipString;
             
            obj.guihandles.offsett=uicontrol('style','text','String','offset (ADU)','Position',[xpos1+2*xw,top-9*vsep,xw*1.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.offset=uicontrol('style','edit','String','100','Position',[xpos1+xw*3.5,top-9*vsep,xw*.5,vsep],'FontSize',fontsize);
            obj.guihandles.offset.TooltipString=sprintf('Offset (ADU).');
            obj.guihandles.offsett.TooltipString=obj.guihandles.offset.TooltipString;
             
            obj.guihandles.peaktext=uicontrol('style','text','String','Peak candidate finding:','Position',[xpos1,top-11*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment',hatitle,'FontWeight','bold');
            
            obj.guihandles.backgroundmodet=uicontrol('style','text','String','Background estimation: ','Position',[xpos1,top-12*vsep,xw*2,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.backgroundmode=uicontrol('style','popupmenu','String',{'Difference of Gaussians (fast)','None'},'Position',[xpos1+xw*2,top-12*vsep,xw*2,vsep],'FontSize',fontsize);
            obj.guihandles.backgroundmode.TooltipString=sprintf('Algorithm for background estimation used for candidate finding.\n Difference of Gaussians is fast');%, wavelet filter is more accurate.');
            obj.guihandles.backgroundmodet.TooltipString=obj.guihandles.backgroundmode.TooltipString;
            
            obj.guihandles.peakfiltert=uicontrol('style','text','String','Filter size (pixel)','Position',[xpos1,top-13*vsep,xw*1.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.peakfilter=uicontrol('style','edit','String','1.2','Position',[xpos1+xw*1.5,top-13*vsep,xw*.5,vsep],'FontSize',fontsize);

            obj.guihandles.peakfilter.TooltipString=sprintf('After background subtraction, the image is filtered with a Gaussain filter prior to peak finding. \n Specify here the width of the Gaussian filter (sigma) in pixels.');
            obj.guihandles.peakfiltert.TooltipString=obj.guihandles.peakfilter.TooltipString;
            
            obj.guihandles.peakcutofft=uicontrol('style','text','String','cutoff (photons)','Position',[xpos1+2*xw,top-13*vsep,xw*1.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.peakcutoff=uicontrol('style','edit','String','1','Position',[xpos1+xw*3.5,top-13*vsep,xw*.5,vsep],'FontSize',fontsize);
            obj.guihandles.peakcutoff.TooltipString=sprintf('Cutoff value to distinguish background from real localizations. \n Units are maximum pixel values of the photon-converted, filtered and background-subtracted image. \n Use the preview function to test several values to find the optimal one.');
            obj.guihandles.peakcutofft.TooltipString=obj.guihandles.peakcutoff.TooltipString;
                      
            obj.guihandles.fittxt=uicontrol('style','text','String','Fitting parameters:','Position',[xpos1,top-15*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment',hatitle,'FontWeight','bold');

            obj.guihandles.roifitt=uicontrol('style','text','String','ROI size (pixel)','Position',[xpos1,top-16*vsep,xw*1.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.roifit=uicontrol('style','edit','String','13','Position',[xpos1+xw*1.5,top-16*vsep,xw*.5,vsep],'FontSize',fontsize);
            obj.guihandles.roifit.TooltipString=sprintf('Size of the fitting region (pixels).'); 
            obj.guihandles.roifitt.TooltipString=obj.guihandles.roifit.TooltipString;
            
            obj.guihandles.bidirectional=uicontrol('style','checkbox','String','2D','Position',[xpos1+2*xw,top-16*vsep,xw*1,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.bidirectional.TooltipString=sprintf('Check this option for a 2D dataset to enable bi-directional fitting.'); 
            obj.guihandles.mirror=uicontrol('style','checkbox','String','mirror','Position',[xpos1+xw*3,top-16*vsep,xw*1,vsep],'FontSize',fontsize);
            obj.guihandles.mirror.TooltipString=sprintf('EMCCD cameras mirror the image if in EM mode. \n Thus, if the bead stack calibration is taken with conventional gain (recommended for better SNR), calibration and experiment do not match. \n In this case, check the mirror option to mirror the single-molcule images prior to fitting.');             
            
            obj.guihandles.outputtxt=uicontrol('style','text','String','Output file:','Position',[xpos1,top-18*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment',hatitle,'FontWeight','bold');
          
            obj.guihandles.selectoutput=uicontrol('style','pushbutton','String','Set output file','Position',[xpos1,top-19*vsep,xw*1.5,vsep],'FontSize',fontsize,'Callback',@obj.selectoutput_callback);
            obj.guihandles.selectoutput.TooltipString='Select output file in which to write the localization table.';
            obj.guihandles.outputfile=uicontrol('style','edit','String','','Position',[xpos1,top-20*vsep,xw*4,vsep],'FontSize',fontsize);
            obj.guihandles.outputfile.TooltipString='Output file that contains the localization table.';
            
            obj.guihandles.outputformat=uicontrol('style','popupmenu','String',{'csv','pointcloud-loader','ViSP'},'Position',[xpos1+1.5*xw,top-19*vsep,xw*2.5,vsep],'FontSize',fontsize);
            obj.guihandles.outputformat.TooltipString=sprintf('Choose output format. CSV saves all fit parameters. Can be opened for instance in PALMsiever. \n Pointcloud-loader and ViSP are popular 3D viewers for localization data.');
                
            obj.guihandles.pixelsizet=uicontrol('style','text','String','pixel size (nm)','Position',[xpos1,top-21*vsep,xw*1.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.pixelsize=uicontrol('style','edit','String','100','Position',[xpos1+xw*1.5,top-21*vsep,xw*.5,vsep],'FontSize',fontsize);
            obj.guihandles.pixelsize.TooltipString='Pixelsize (nm). Required for some output formats that require all units in nm.';
            obj.guihandles.pixelsizet.TooltipString=obj.guihandles.pixelsize.TooltipString;
            
            obj.guihandles.preview=uicontrol('style','pushbutton','String','Preview frame:','Position',[xpos1,top-23*vsep,xw*1.5,vsep],'FontSize',fontsize, 'Callback',@obj.preview_callback);
            obj.guihandles.preview.TooltipString='Localize the frame with the number specified here to check all fitting parameters.';
            
            obj.guihandles.previewframe=uicontrol('style','edit','String','1','Position',[xpos1+1.5*xw,top-23*vsep,xw*.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.previewframe.TooltipString='Frame used for preview. It is recommended to check several frames from different parts of the data set.';
            
            obj.guihandles.localize=uicontrol('style','pushbutton','String','Localize','Position',[xpos1+2.5*xw,top-23*vsep,xw*1.5,vsep],'FontSize',fontsize, 'Callback',@obj.localize_callback,'FontWeight','bold');
            obj.guihandles.localize.TooltipString='Perform 3D fitting of single molecules';

            obj.guihandles.status=uicontrol('style','text','String','Status','Position',[xpos1,top-25*vsep,xw*3.5,vsep],'FontSize',fontsize);
            
            obj.guihandles.stop=uicontrol('style','togglebutton','String','Stop','Position',[xpos1+3.5*xw,top-25*vsep,xw*0.5,vsep],'FontSize',fontsize, 'Callback',@obj.stop_callback);
            obj.guihandles.stop.TooltipString='Perform 3D fitting of single molecules';
        end
        function selectfiles_callback(obj,a,b,which)
            switch which
                case 1
                    switch obj.guihandles.loader.Value
                        case {1,2}
                            file=obj.guihandles.imagefile.String;
                            if isempty(file)
                                file='*.tif';
                            end
                            handle='imagefile';
                        case 3 %ImageJ
                            fijipath='ImageJ';
                            if ~isempty(obj.mij)
                                try
                                obj.mij.exit;
                                catch err
                                    err
                                end
                            end
                            wf=msgbox('Open the raw image file in ImageJ. Then set parameters, load 3D calibratin and use preview/localize. Don''t close ImageJ.');
                            waitfor(wf);
                            myMiji(true,fijipath);
                            obj.mij=MIJ;
                             obj.guihandles.imagefile.String='from_ImageJ';
                            return
                    end
                case 2
                    file=obj.guihandles.calfile.String;
                    if isempty(file)
                        file='*.mat';
                    end
                    handle='calfile';
                 case 3
                    file=obj.guihandles.scmosfile.String;
                    if isempty(file)
                        file='*.*';
                    end
                    handle='scmosfile';                   
            end
            [f,p]=uigetfile(file);
            if f
                obj.guihandles.(handle).String=[p f];
                if which==1
                    if isempty(obj.guihandles.outputfile.String)
                        obj.guihandles.outputfile.String=strrep([p f],'.tif','.csv');
                    end
                end
            end
 
        end
        function selectoutput_callback(obj,a,b)
            of=obj.guihandles.outputfile.String;
            if isempty(of)
           
                of='out.csv';
            end
             if contains(obj.guihandles.outputformat.String{obj.guihandles.outputformat.Value},'ViSP')
                 of=strrep(of,'.csv','.3dlp');
             else
                 of=strrep(of,'.3dlp','.csv');
             end
            [f,p]=uiputfile(of);
            if f
            obj.guihandles.outputfile.String=[p,f];
            end
        end
        
        
        function preview_callback(obj,a,b)
            obj.guihandles.stop.Value=false;
            obj.stop_callback(obj.guihandles.stop);
            p=obj.getguiparameters;
            p.status.String='Preview...'; drawnow;
            p.preview=true;
            simplefitter_cspline(p)
        end
        
        
        function localize_callback(obj,a,b)
            obj.guihandles.stop.Value=false;
            obj.stop_callback(obj.guihandles.stop);
            p=obj.getguiparameters;
            p.preview=false;
            
            if isempty(p.outputfile)
                errordlg('please define output file');
                return
            end
            p.status.String='Start localization...';drawnow;
            simplefitter_cspline(p)
        end        
        
        function p=getguiparameters(obj)
            %read the parameters of the GUI and write them into parameter
            %structure
            
            p.imagefile=obj.guihandles.imagefile.String;
            p.calfile=obj.guihandles.calfile.String;
            p.offset=str2double(obj.guihandles.offset.String);
            p.conversion=str2double(obj.guihandles.conversion.String);
            p.previewframe=str2double(obj.guihandles.previewframe.String);
            p.peakfilter=str2double(obj.guihandles.peakfilter.String);
            p.peakcutoff=str2double(obj.guihandles.peakcutoff.String);
            p.roifit=str2double(obj.guihandles.roifit.String);
            p.bidirectional=(obj.guihandles.bidirectional.Value);
            p.mirror=(obj.guihandles.mirror.Value);
            p.status=obj.guihandles.status;
            p.outputfile=obj.guihandles.outputfile.String;
            p.outputformat=obj.guihandles.outputformat.String{obj.guihandles.outputformat.Value};
            p.pixelsize=str2num(obj.guihandles.pixelsize.String);
            p.loader=(obj.guihandles.loader.Value);
            p.mij=obj.mij;
            p.backgroundmode=obj.guihandles.backgroundmode.String{obj.guihandles.backgroundmode.Value};
            
            
            p.isscmos=obj.guihandles.isscmos.Value;
            p.scmosfile=obj.guihandles.scmosfile.String;
            
            if p.isscmos && p.mirror
                warndlg('sCMOS cameras do not mirror image, this only happens for EMCCD cameras when bead calibration is taken without and image data are taken with EM gain. Mirror not used.')
                p.mirror=false;
            end
            
%             if isempty(p.calfile)
%                 warndlg('3D bead calibration file not selected. Using Gaussian fitter');
%                 p.isspline=false;
%             else
%                 p.isspline=true;
% %                 errordlg('please load 3D bead calibration file first')
% %                 error('please load 3D bead calibration file first')
%             end
            if isempty(p.imagefile)
                errordlg('please load image file first')
                error('please load image file first')
            end
        end
        
        function changeloader_callback(obj,a,b)
        end
        
        function stop_callback(obj,object,b)
            global simplefitter_stop
            simplefitter_stop=object.Value;
        end
       
    end
end
