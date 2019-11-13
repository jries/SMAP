classdef PeakFinder<interfaces.WorkflowModule
%     Performs peak-finding according either to local maximum finding  or
%     to non-maximum suppression. Cutoff based on absolute numbers,
%     probabilities, or dynamically calculated difference from background.
%     (A. Neubeck and L. Van Gool, ?Efficient non-maximum suppression,?
%     presented at the 18th International Conference on Pattern
%     Recognition, Vol 3, Proceedings, 10662 LOS VAQUEROS CIRCLE, PO BOX
%     3014, LOS ALAMITOS, CA 90720-1264 USA, 2006, pp. 850?855.).';
% 
    properties (Access=private)
        probability=0.05;
        dynamicfactor=1.7;
        absolutecutoff=1;
        roimask
        preview
    end
    methods
        function obj=PeakFinder(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
            obj.inputParameters={'loc_loc_filter_sigma','EMon'};
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.guihandles.cutoffmode.Callback={@cutoffmode_callback,obj};
            obj.guihandles.cutoffvalue.Callback={@cutoffvalue_callback,obj};
             obj.guihandles.peakfindmethod.Callback={@peakfindmethod_callback,obj};
             cutoffvalue_callback(0,0,obj)
        end
        function prerun(obj,p)
            cutoffvalue_callback(0,0,obj)
            obj.roimask=obj.getPar('loc_roimask');
             obj.preview=obj.getPar('loc_preview');
        end
        function dato=run(obj,data,p)
            image=data.data;%get;
            if ~isempty(image)
                if all(size(obj.roimask)== size(image))
                    image(~obj.roimask)=-1;
                end
%             switch p.peakfindmethod.Value
%                 case 1 %maximum
                    maxima=maximumfindcall(image); %find maxima
                    maxima(:,[1 2])=maxima(:,[2 1]);
                    dynamicf=1;
%                 case 2 %NMS
%                 maxima=NMS2DBlockCcall(image,p.NMS_kernel_size); %find maxima
%                 dynamicf=3/p.NMS_kernel_size;
%             end
                switch p.cutoffmode.Value
                    case 1
                        cutoff=getdynamiccutoff(maxima,obj.dynamicfactor*dynamicf);
                    case {2,3} 
                        cutoff=obj.absolutecutoff;
                end
    %             cutoff
                maxind= (maxima(:,3)>cutoff);

                if p.use_mindistance
                    maxima=maxima(maxind,:);
                    maxind=~tooclose(maxima(:,1),maxima(:,2),p.mindistance); % p.mindistance) & ~tooclose(,p.mindistance);
                end
                maxout.ypix=maxima(maxind,1);
                maxout.xpix=maxima(maxind,2);
                maxout.phot=maxima(maxind,3);
                maxout.frame=0*maxout.xpix+data.frame;
                dato=data;         
                dato.data=(maxout); 
                if obj.preview
                    obj.setPar('preview_peakfind',maxout);
                end
            else

                dato=data;
                dato.data=[];
            end
        end
    end
end

function val = prob2photon(p,PSFx0,filtersize,excess)
val=sqrt(2)*erfinv(1-2*p); %in paper it is 2p-1 
 val=val*sqrt(PSFx0^2/(PSFx0^2+filtersize^2)*excess);
 val=max(1E-7,val);
val=min(1E7,val);
if isempty(val)
    val=1;
end
end

function p=photon2prob(val,PSFx0,filtersize,excess)
val=val/sqrt(PSFx0^2/(PSFx0^2+filtersize^2)*excess);
p=(1-erf(val/sqrt(2)))/2;
p=max(1E-7,p);
p=min(1,p);
if isempty(p)
    p=0.01;
end
end

% function peakfindmethod_callback(a,b,obj)
% if obj.guihandles.peakfindmethod.Value==1
%     obj.guihandles.NMS_kernel_size.Visible='off';
% else
%     obj.guihandles.NMS_kernel_size.Visible='on';
% end
% end

function cutoffmode_callback(a,b,obj)
p=obj.getAllParameters;
state=p.cutoffmode.Value; 
 switch state    
     case 1 %dynamic
         pset.cutoffvalue=obj.dynamicfactor;
     case 2 %probability
         pset.cutoffvalue=obj.probability;
     case 3 %absolute
         pset.cutoffvalue=obj.absolutecutoff;
end
obj.setGuiParameters(pset);
end

function cutoffvalue_callback(a,b,obj)
p=obj.getAllParameters;
PSFx0=1;
state=p.cutoffmode.Value;
excess=double(p.EMon)+1;
if isempty(excess)
    excess=1;
end
switch state
    case 1 %dynamic
        obj.dynamicfactor=p.cutoffvalue;
    case 2 %p       
        obj.probability=p.cutoffvalue;
        obj.absolutecutoff=prob2photon(p.cutoffvalue,PSFx0,p.loc_loc_filter_sigma(1),excess);
    case 3 %abs      
        obj.absolutecutoff=p.cutoffvalue;
        obj.probability=photon2prob(p.cutoffvalue,PSFx0,p.loc_loc_filter_sigma(1),excess);
end
end

function co=getdynamiccutoff(maxima,factor)
ps=[.2 .5 .8];
if size(maxima,1)<10
    if isempty(maxima)
        co=0;
    else
        co=mean(maxima(:,3))*factor;
    end
else
    
qs=myquantilefast(maxima(:,3),ps);
slope=(qs(3)-qs(1))/(ps(3)-ps(1));
co=qs(2)+slope*0.5*2*factor;
end
end

function pard=guidef
pard.cutoffstring.object=struct('Style','text','String','       Cutoff: ');
pard.cutoffstring.position=[1,1];

pard.cutoffmode.object=struct('Style','popupmenu','String',{{'dynamic (factor)','probability (p<1)','absolute (photons)'}},'Value',1);
pard.cutoffmode.position=[1,1];
pard.cutoffmode.Width=1.5;
pard.cutoffmode.TooltipString=sprintf('How to determine the cutoff: \n Dynamic: use the distribution of pixel intensity to estimate likely localizations. Factor: adjust sensitivity. \n Probability: use probabilistic model (SimpleSTORM) to determine the likelyhood for pixel being localization. \n Absolut: Pixel intensity in normalized image. \n Choose display=normalized to read out thes normalized values.');
pard.cutoffmode.Optional=true;

pard.cutoffvalue.object=struct('Style','edit','String','1.7');
pard.cutoffvalue.position=[1,2.5];
pard.cutoffvalue.Width=.5;
pard.cutoffvalue.TooltipString=sprintf('Dynamic: relative factor. \n Probability: directly probability p. \n absolute: cutoff in photons');

pard.use_mindistance.object=struct('Style','checkbox','String','minimum distance (pix)');
pard.use_mindistance.position=[2,1];
pard.use_mindistance.Width=1.5;
pard.mindistance.object=struct('Style','edit','String','7','Visible','on');
pard.mindistance.position=[2,2.5];
pard.mindistance.Width=.5;
% pard.peakfindmethod.object=struct('Style','popupmenu','String',{{'maximum','NMS: kernel size (pix)'}});
% pard.peakfindmethod.position=[2,1];
% pard.peakfindmethod.Width=1.5;
% pard.peakfindmethod.TooltipString=sprintf('Maximum: all local maxima. \n NMS: non-maxiumum suprression. Finds only maxima spaced at least NMS size.');
% pard.peakfindmethod.Optional=true;
% pard.NMS_kernel_size.object=struct('Style','edit','String','5','Visible','off');
% pard.NMS_kernel_size.position=[2,2.5];
% pard.NMS_kernel_size.Width=.5;
% pard.NMS_kernel_size.Optional=true;
pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Performs peak-finding according either to local maximum finding  or to non-maximum suppression. Cutoff based on absolute numbers, probabilities, or dynamically calculated difference from background. (A. Neubeck and L. Van Gool, ?Efficient non-maximum suppression,? presented at the 18th International Conference on Pattern Recognition, Vol 3, Proceedings, 10662 LOS VAQUEROS CIRCLE, PO BOX 3014, LOS ALAMITOS, CA 90720-1264 USA, 2006, pp. 850?855.).';
end