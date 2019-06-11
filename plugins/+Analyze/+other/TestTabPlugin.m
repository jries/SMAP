classdef TestTabPlugin<interfaces.DialogProcessor
% study how to use tabs in module interfaces. 
    properties
        
    end
    methods
        function obj=TestTabPlugin(varargin)    
            obj@interfaces.DialogProcessor(varargin{:}) ;
        end
%         function makeGui(obj,varargin)
%             Vrimold=obj.guiPar.Vrim;handleold=obj.handle;
%            
%             if nargin >2 && varargin{2}==true
%                  obj.guiPar.Vrim=50;
%                  makeGui@interfaces.GuiModuleInterface(obj,varargin{1});
%             else
%                 makeGui@interfaces.DialogProcessor(obj,varargin{:});
%                 obj.guiPar.Vrim=50;
%                 obj.guihandles.tabgroup=uitabgroup(obj.handle,'Position',[0 0 1 .75],'SelectionChangedFcn',{@selectLayer_callback,obj});
%                 obj.guihandles.tab1=uitab(obj.guihandles.tabgroup,'Title','Main');
%                 obj.handle=obj.guihandles.tab1;
%                 makeGui@interfaces.GuiModuleInterface(obj,guidef1);
%                 
%                 obj.addguitotab(1)
%                 
%                 obj.guihandles.tab3=uitab(obj.guihandles.tabgroup,'Title','+');
%             end
%             obj.handle=handleold;
%             obj.guiPar.Vrim=Vrimold;
%         end
        function out=run(obj,p)  
           p
        end
        function pard=guidef(obj)
            pard=guidef1;
            pard.plugininfo.name='Test Tab Gui';
            pard.plugininfo.description= 'x';
            pard.plugininfo.type='ProcessorPlugin';
        end
        function addguitotab(obj,number)
            tag=['Module' num2str(number)];
            obj.guihandles.(['tab' num2str(number)])=uitab(obj.guihandles.tabgroup,'Title',tag,'Tag',tag);
            guidefhere=addnumbertofield(guidefmodule,number);
            Vrimold=obj.guiPar.Vrim;handleold=obj.handle;
            obj.guiPar.Vrim=50;
            obj.handle=obj.guihandles.(['tab' num2str(number)]);
            obj.makeGui(guidefhere,1);
            obj.handle=handleold;
            obj.guiPar.Vrim=Vrimold;
        end
    end
end



function selectLayer_callback(tabgroup,eventdata,obj)
layertitle=(eventdata.NewValue.Title);
if strcmp(layertitle,'+')
    numtabs=length(tabgroup.Children);
    obj.addguitotab(numtabs-1)
    s=1:length(tabgroup.Children);
    s(end-1)=s(end);
    s(end)=s(end)-1;
    tabgroup.Children=tabgroup.Children(s);
    tabgroup.SelectedTab=tabgroup.Children(end-1); 

end
end

function out=addnumbertofield(in,number)
fn=fieldnames(in);
for k=1:length(fn)
    out.([fn{k} '_' num2str(number)])=in.(fn{k});
end
end

function pard=guidef1

pard.tab.tab1='first';
pard.tab.tab2='2';

pard.ch1t.object=struct('String','Ch 1','Style','text');
pard.ch1t.position=[2,1];
pard.ch1t.Width=0.3;
pard.ch1.object=struct('String',{{'layer 1','layer 2','layer 3'}},'Style','popupmenu');
pard.ch1.position=[2,1.3];
pard.ch1.Width=1;
pard.ch1.tab='tab2';

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


function pard=guidefmodule
pard.ch1.object=struct('String','text','Style','text');
pard.ch1.position=[2,1];
pard.ch1.Width=0.3;

pard.rrange.object=struct('String','50','Style','edit');
pard.rrange.position=[3,2];
pard.rrange.Width=.5;
end