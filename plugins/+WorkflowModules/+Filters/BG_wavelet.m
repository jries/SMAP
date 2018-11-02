classdef BG_wavelet<interfaces.WorkflowModule
    properties
        gpufit;
        cutoff;
    end
    methods
        function obj=BG_wavelet(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
            obj.outputChannels=1; %1: image-background. 2: background image
            obj.outputParameters={'loc_subtractbg'};
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.setPar('loc_blocksize_frames',1); %preview: tiff loader should load only one frame
            filtermode_callback(obj.guihandles.filtermode,0,obj)
        end
        function prerun(obj,p)
            if p.filtermode.Value>2  
                try
                    gpuArray(5);
                    obj.gpufit=true;
                catch err
                    disp('cuda did not work')
                    err
                    obj.gpufit=false;
                end
            else
                obj.gpufit=false;
            end
        end
        function output=run(obj,data,p)
           
            if ~p.loc_subtractbg
                dato2=data;
                dato2.data=(0*data.data);
                output=dato2;
                return
            end
            if ~isempty(data.data)
                img=data.data;
                
                if true %pre-filter
                    co=myquantilefast(img(:),.99);
                    img(img>co)=co;
                end
                if p.filtermode.Value>2&&obj.gpufit
                    img=gpuArray(img);
                end
                switch p.filtermode.Value
                    case 1 %wavelet
                        bg=mywaveletfilter(img,p.loc_wavelet_level,p.loc_wavelet_refine,true);
%                          bg=mywaveletfilter(gpuArray(img),p.loc_wavelet_level,p.loc_wavelet_refine,true);
                    case 2 %minimum
                      
                        bg=minfilt2_min_fast(img,p.min_filter_size);
                    case 3 % a trous
%                         imgg=gpuArray((img));
                        bg=(mywaveletfilteratrous(img,false));

%                         bg=gather(bgg);
%                         clear imgg;
                end
                if p.filtermode.Value>2&&obj.gpufit
                    bg=gather(bg);
                    clear img;
                end
                dato=data;
                dato.data=bg;
                output=dato; 
            else %eof
                output=data;
            end 
        end
       
    end
end

function filtermode_callback(object,b,obj)

switch object.Value
    case 1
        onf={'loc_wavelet_levelt','loc_wavelet_level','loc_wavelet_refine'};
        offf={'min_filter_sizet','min_filter_size'};
    case 2
        offf={'loc_wavelet_levelt','loc_wavelet_level','loc_wavelet_refine'};
        onf={'min_filter_sizet','min_filter_size'};
    case 3
        offf={'min_filter_sizet','min_filter_size','loc_wavelet_levelt','loc_wavelet_level','loc_wavelet_refine'};
        onf={};
end

obj.fieldvisibility('on',onf,'off',offf);

end

function pard=guidef(obj)
pard.loc_subtractbg.object=struct('Style','checkbox','String','Subtract background','Value',1);
pard.loc_subtractbg.position=[1,1];
pard.loc_subtractbg.Width=2;
pard.loc_subtractbg.TooltipString=sprintf('If checked, the background is subtracted for Peak finding. \n This does NOT mean, that fitting is performed on the background corrected images.');
% pard.text1.object=struct('Style','text','String','Wavelet filtering:');
% pard.text1.position=[2,1];
% pard.text1.Optional=true;


pard.filtermode.object=struct('Style','popupmenu','String',{{'Wavelet','fast minimum filter', 'a trous'}},'Value',1,'Callback',{{@filtermode_callback,obj}});
pard.filtermode.position=[2,1];
pard.filtermode.Width=2;
pard.filtermode.TooltipString=sprintf('trous algorithm is used. Only level 2. Can be faster on GPU.');
pard.filtermode.Optional=true;

% pard.loc_wavelet_atrous.object=struct('Style','checkbox','String','fast a trous','Value',0,'Callback',{{@atrous_callback,obj}});
% pard.loc_wavelet_atrous.position=[2,2];
% pard.loc_wavelet_atrous.Width=1.3;
% pard.loc_wavelet_atrous.TooltipString=sprintf('If checked, the a trous algorithm is used. Only level 2. Can be faster on GPU.');
% pard.loc_wavelet_atrous.Optional=true;

pard.min_filter_sizet.object=struct('Style','text','String','filter size');
pard.min_filter_sizet.position=[3,1.3];
pard.min_filter_sizet.Optional=true;
pard.min_filter_size.object=struct('Style','edit','String','5');
pard.min_filter_size.position=[3,2.3];
pard.min_filter_size.Width=.7;
pard.min_filter_size.TooltipString=sprintf('Size of kernel of minimum filter in pixels');
pard.min_filter_size.Optional=true;


pard.loc_wavelet_levelt.object=struct('Style','text','String','Wavelet level');
pard.loc_wavelet_levelt.position=[3,1.3];
pard.loc_wavelet_levelt.Optional=true;


pard.loc_wavelet_level.object=struct('Style','edit','String','3');
pard.loc_wavelet_level.position=[3,2.3];
pard.loc_wavelet_level.Width=.7;
pard.loc_wavelet_level.TooltipString=sprintf('Wavelet level for background correction. Typical: 3 (range: 2-6)');
pard.loc_wavelet_level.Optional=true;


pard.loc_wavelet_refine.object=struct('Style','checkbox','String','Refined background estimation','Value',0);
pard.loc_wavelet_refine.position=[4,1.3];
pard.loc_wavelet_refine.Width=2;
pard.loc_wavelet_refine.TooltipString=sprintf('Iterative refinement of background estimation. \n Slower, use mainly for: \n a) very bright fluorophores (e.g. beads) that otherwise lead to ghost localizatiosn, \n b) for high background and detection of weak fluorophores.');
pard.loc_wavelet_refine.Optional=true;


pard.plugininfo.type='WorkflowModule';
pard.plugininfo.description='Wavelet based background estimation. Two implementeations: A trous (can be faster on GPU, but only level 2, adepted from I. Izeddin, J. Boulanger, V. Racine, C. G. Specht, A. Kechkar, D. Nair, A. Triller, D. Choquet, M. Dahan, and J. B. Sibarita, ?Wavelet analysis for single molecule localization microscopy,? Opt Express, vol. 20, no. 3, pp. 2081?2095, Jan. 2012.) and a direct discreet wavelet transform.';
end