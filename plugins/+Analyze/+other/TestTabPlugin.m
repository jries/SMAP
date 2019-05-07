classdef TestTabPlugin<interfaces.DialogProcessor
%  pair correlation functions calculated according to:
%  Sengupta, Prabuddha, Tijana Jovanovic-Talisman, Dunja Skoko, Malte Renz, 
%  Sarah L Veatch, and Jennifer Lippincott-Schwartz. “Probing Protein Heterogeneity 
%  in the Plasma Membrane Using PALM and Pair Correlation Analysis.” 
%  Nature Methods 8 (September 18, 2011): 969.
    properties
        
    end
    methods
        function obj=TestTabPlugin(varargin)    
            obj@interfaces.DialogProcessor(varargin{:}) ;
        end
        function makeGui(obj,varargin)
             Vrimold=obj.guiPar.Vrim;
             makeGui@interfaces.DialogProcessor(obj,varargin{:});
    
            obj.guiPar.Vrim=50;
            obj.guihandles.tg=uitabgroup(obj.handle,'Position',[0 0 1 .75],'SelectionChangedFcn',{@selectLayer_callback,obj});
            handleold=obj.handle;
            obj.guihandles.tab1=uitab(obj.guihandles.tg,'Title','Main');
            obj.handle=obj.guihandles.tab1;
            makeGui@interfaces.GuiModuleInterface(obj,guidef1);
            obj.guihandles.tab2=uitab(obj.guihandles.tg,'Title','Module 1');
            obj.handle=obj.guihandles.tab2;
            makeGui@interfaces.GuiModuleInterface(obj,guidef2);
            obj.guihandles.tab3=uitab(obj.guihandles.tg,'Title','+');
            obj.handle=handleold;
            obj.guiPar.Vrim=Vrimold;
        end
        function out=run(obj,p)  
           p
        end
        function pard=guidef(obj)
            pard.plugininfo.name='Test Tab Gui';
            pard.plugininfo.description= 'x';
            pard.plugininfo.type='ProcessorPlugin';
        end
%         function pard=guidef2(obj)
%             pard=guidef2;
%         end

    end
end


function selectLayer_callback(tabgroup,eventdata,obj)
layer=(eventdata.NewValue.Tag);
layertitle=(eventdata.NewValue.Title);
if strcmp(layertitle,'+')
    newlayernumber=obj.numberOfLayers+1;
    tag=['Layer' num2str(newlayernumber)];
      h.(['tab_layer' num2str(newlayernumber)])=uitab(obj.guihandles.layertab,'Title',tag,'Tag',tag);
       h.(['layer_' num2str(newlayernumber)])=obj.addlayer(  h.(['tab_layer' num2str(newlayernumber)]),newlayernumber);
    obj.numberOfLayers=newlayernumber;
    obj.setPar('numberOfLayers',newlayernumber);
    s=1:length(tabgroup.Children);
    s(end-1)=s(end);
    s(end)=s(end)-1;
    tabgroup.Children=tabgroup.Children(s);
    tabgroup.SelectedTab=tabgroup.Children(end-1); 
    selectedlayer=newlayernumber;
    obj.setPar('layercheck',true,'Value','layer',selectedlayer)
else
    selectedlayer=sscanf(layer,'Layer%d');
end
for k=1:length(tabgroup.Children)-1
%     layernumber=tabgroup.Children(k).Title;
    layertag=tabgroup.Children(k).Tag;
    obj.children.(layertag).setlayer(selectedlayer);
end
updatelayernames(obj)
end


function pard=guidef1
pard.ch1t.object=struct('String','Ch 1','Style','text');
pard.ch1t.position=[2,1];
pard.ch1t.Width=0.3;
pard.ch1.object=struct('String',{{'layer 1','layer 2','layer 3'}},'Style','popupmenu');
pard.ch1.position=[2,1.3];
pard.ch1.Width=1;
pard.ch2t.object=struct('String','Ch 2','Style','text');
pard.ch2t.position=[2,3];
pard.ch2t.Width=0.3;
pard.ch2.object=struct('String',{{'none','layer 1','layer 2','layer 3'}},'Style','popupmenu');
pard.ch2.position=[2,3.3];
pard.ch2.Width=1;

pard.rranget.object=struct('String','Range r (nm)','Style','text');
pard.rranget.position=[3,1];
pard.rranget.Width=1;
pard.rrange.object=struct('String','1000','Style','edit');
pard.rrange.position=[3,2];
pard.rrange.Width=.5;

pard.rpixt.object=struct('String','Pixelsize (nm)','Style','text');
pard.rpixt.position=[3,3];
pard.rpixt.Width=1;
pard.rpix.object=struct('String','2','Style','edit');
pard.rpix.position=[3,4];
pard.rpix.Width=.5;



end


function pard=guidef2
pard.ch1t22.object=struct('String','text','Style','text');
pard.ch1t22.position=[2,1];
pard.ch1t22.Width=0.3;

pard.rrange22.object=struct('String','50','Style','edit');
pard.rrange22.position=[3,2];
pard.rrange22.Width=.5;




end