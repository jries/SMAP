classdef MedianBGcalculator<interfaces.WorkflowModule
    properties
        blockOfFrames=single(0);
        blockIndex;
%         datablock=interfaces.WorkflowData;
        runparameters
    end
    methods
        function obj=MedianBGcalculator(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
            obj.outputChannels=2; %1: image-background. 2: background image
            obj.outputParameters={'loc_bg_dx','loc_blocksize_frames','loc_subtractbg'};
        end
        function pard=guidef(obj)
            pard=guidef;
        end
%         function initGui(obj)
%             initGui@interfaces.WorkflowModule(obj);
%         end
        function prerun(obj,p)
            obj.blockIndex=[];
           obj.runparameters=p;
           if obj.runparameters.loc_interleavedexc
               obj.runparameters.numberOfPulses=2;
           else
               obj.runparameters.numberOfPulses=1;
           end
        end
        function output=run(obj,data,p)
            persistent block datablock 
            output=[];
            par=obj.runparameters;
            if ~par.loc_subtractbg
                obj.output(data,2)
                dato2=data;
                dato2.data=(0*data.data);
                obj.output(dato2,1);
%                 output=dato2;
                return
            end
            img=data.data;
            blockIndexh=obj.blockIndex;
            dt=par.loc_blocksize_frames;
            dx=par.loc_bg_dx;

            %initialization
            if isempty(blockIndexh) %first time
                s=size(img);
                for k=1:par.numberOfPulses
                    block{k}=zeros(s(1),s(2),dt,'single');
                    datablock{k}(dt)=data;%.copy;
                end
                blockIndexh=zeros(par.numberOfPulses,par.numberOfPulses);

            end
                
            if ~data.eof
                frame=data.frame;
                pulsenumber=mod(frame,par.numberOfPulses)+1;

                %fill in block
                blockIndexh(pulsenumber)=blockIndexh(pulsenumber)+1;
                blockIndexP=blockIndexh(pulsenumber);
%                 ind=blockIndex(pulsenumber);
                block{pulsenumber}(:,:,blockIndexP)=img;
                datablock{pulsenumber}(blockIndexP)=data;

                 %do calculatations
%                  blockIndexP
                if blockIndexP==dt
                    if dx<3
                        bgIm=mymedian(block{pulsenumber}(:,:,1:blockIndexP),3); %only temporal median along dimension 3   
                    else
                        bgIm=mymedianfilter(block{pulsenumber}(:,:,1:blockIndexP),dx);
                    end
                    bgIm(isnan(bgIm))=0;
                    for k=1:blockIndexP
                        try
                        dato=datablock{pulsenumber}(k);

                        dato.data=bgIm;  
                        obj.output(dato);
                        catch
                            disp('problem')
                        end
                    end
                    blockIndexh(pulsenumber)=0;
                end
            else %eof
                %output the rest of the data
                for k=1:par.numberOfPulses
                    if blockIndexh(k)>0
                        if dx<3
                            bgIm=mymedian(block{k}(:,:,1:blockIndexh(k)),3); %only temporal median along dimension 3   
                        else
                            bgIm=mymedianfilter(block{k}(:,:,1:blockIndexh(k)),dx);
                        end
                        bgIm(isnan(bgIm))=0;
                        for l=1:blockIndexh(k)
                            dato=datablock{k}(l);
                            dato.data=bgIm;  
                            obj.output(dato);
                        end
                        blockIndexh(k)=0;
                    end
                end
                obj.output(data,1);
            end
            obj.blockIndex=blockIndexh;
        end       
    end
end


function pard=guidef
pard.loc_subtractbg.object=struct('Style','checkbox','String','Subtract background','Value',1);
pard.loc_subtractbg.position=[1,1];
pard.loc_subtractbg.Width=2;
pard.loc_subtractbg.TooltipString=sprintf('If checked, the background is subtracted for Peak finding. \n This does NOT mean, that fitting is performed on the background corrected images.');

pard.loc_interleavedexc.object=struct('Style','checkbox','String','Interleaved ','Value',0);
pard.loc_interleavedexc.position=[2,2];
pard.loc_interleavedexc.TooltipString=sprintf('Check for interleaved excitation with two lasers to calculate the BG only with corresponding frames.');

pard.text1.object=struct('Style','text','String','Median filtering:');
pard.text1.position=[2,1];

pard.text2.object=struct('Style','text','String','dx (pixels)');
pard.text2.position=[3,1.3];
pard.loc_bg_dx.object=struct('Style','edit','String','3');
pard.loc_bg_dx.position=[3,2.3];
pard.loc_bg_dx.Width=.7;
pard.loc_bg_dx.TooltipString=sprintf('Number of pixels over which the median is calculated. (0-5)');

pard.text3.object=struct('Style','text','String','dt (frames)');
pard.text3.position=[4,1.3];


pard.loc_blocksize_frames.object=struct('Style','edit','String','100');
pard.loc_blocksize_frames.position=[4,2.3];
pard.loc_blocksize_frames.Width=.7;
pard.loc_blocksize_frames.TooltipString=sprintf('Number of frames over which the median is calculated. (Typically 100. range: 0-1000)');
pard.plugininfo.type='WorkflowModule';
pard.plugininfo.description='Background calcualtion based on median filtering. For reasonable computational complexity, the median is calculated for blocks of dx*dx pixels and dt frames. Values are interpolated in x and y.';
end