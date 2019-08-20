classdef modelFitter<interfaces.SEEvaluationProcessor
%     Fit parametrized model to SMLM data. Not public yet.
    properties

    end
    methods
        function obj=modelFitter(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, inp)
            % Import the kimograph
%             site = obj(1).locData
            layers=find(inp.sr_layerson);
            locs=obj.getLocs({'locprecnm','PSFxnm','xnm','ynm','frame'},'layer',1,'size',inp.se_siteroi);      
        end
        function makeGui(obj,varargin)
            if nargin >2 && varargin{2}==true %gui of tabs
                 makeGui@interfaces.GuiModuleInterface(obj,varargin{1});
            else %real make GUI of main gui
                makeGui@interfaces.GuiModuleInterface(obj); %make the main GUI from the guidef definitions
                %Settings               
                optimizernames={'f1','f2'};%if possible, get from fitter object
                obj.guihandles.optimizer.String=optimizernames;
                Data={'Name',15};
                oldh=obj.guihandles.optimizerpar;
                pos=oldh.Position;
                htable=uitable(oldh.Parent,'Data',Data,'Position',pos);
                delete(oldh);
                obj.guihandles.optimizerpar=htable;

                %add first module and  + tab
                obj.guihandles.tabgroup=obj.guihandles.tab1.Parent;
                obj.addguitotab(1);
                obj.guihandles.addtab=uitab(obj.guihandles.tabgroup,'Title','+');
                obj.guihandles.tabgroup.SelectionChangedFcn={@selectLayer_callback,obj};
                
            end
            
            
        end
     
        function pard=guidef(obj)
            pard=guidef;
        end
        function addguitotab(obj,number)
            %this adds model (number) to the tab group
            tag=['M' num2str(number)];
            obj.guihandles.(['tab' num2str(number)])=uitab(obj.guihandles.tabgroup,'Title',tag,'Tag',tag);
            guidefhere=addnumbertofield(guidefmodel(obj,number),number);
            Vrimold=obj.guiPar.Vrim;handleold=obj.handle;
            obj.guiPar.Vrim=0;
            obj.handle=obj.guihandles.(['tab' num2str(number)]);
            obj.makeGui(guidefhere,1);
            obj.handle=handleold;
            obj.guiPar.Vrim=Vrimold;
            
            %initialize conversion table. Here change the looks!
            hconv=obj.guihandles.(['tabconv_' num2str(number)]);
            dataconv={'field1','expression','field2'};
            obj.guihandles.(['convtable_' num2str(number)])=uitable(hconv,'Data',dataconv);
            obj.guihandles.(['convtable_' num2str(number)]).Position(3:4)=obj.guihandles.(['convtable_' num2str(number)]).Position(3:4)-[50 150];
            
            %initialize parameter table. Here change the looks!
            hpar=obj.guihandles.(['tabpar_' num2str(number)]);
            datapar={'Name','fix','Value','LB','UB'};
            obj.guihandles.(['partable_' num2str(number)])=uitable(hpar,'Data',datapar);
            obj.guihandles.(['partable_' num2str(number)]).Position(3:4)=obj.guihandles.(['partable_' num2str(number)]).Position(3:4)-[50 150];
        end
    end

end

function selectLayer_callback(tabgroup,eventdata,obj) 
% if + tab selected this makes a new model
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
%helper function. obj.guihandles needs to be flat. Add number to field
%name.
fn=fieldnames(in);
for k=1:length(fn)
    if strcmp(fn{k},'tab')
        out.tab=addnumbertofield(in.(fn{k}),number);
    else
    out.([fn{k} '_' num2str(number)])=in.(fn{k});
    end
end
end

function loadmodel_callback(a,b,obj)
% executed when model is loaded
modelnumber=(a.Parent.Parent.Parent.Title(2:end)); %hack to get the right tab
fnold=obj.guihandles.(['modelname_' modelnumber]).String;
[f,p]=uigetfile(fnold);
if ~f %no model selected: return
    return
end
obj.guihandles.(['modelname_' modelnumber]).String=[p f];
initmodel(obj, modelnumber)
end

function initmodel(obj, modelnumber)
%after the model is selected, you need to update the model parameter table.
htab=obj.guihandles.(['partable_' modelnumber]); %handle of table
htab.Data;
end

function pard=guidef(obj)
pard.tab.tab1='Settings';
pard.optimizer.object=struct('Style','popupmenu','String',{{'text','t2'}});
pard.optimizer.position=[2,1];
pard.optimizer.Width=3.5;
pard.optimizer.tab='tab1';
pard.optimizerpar.object=struct('Style','text','String','');
pard.optimizerpar.position=[6,1];
pard.optimizerpar.Width=3.5;
pard.optimizerpar.Height=4;
pard.optimizerpar.tab='tab1';

pard.loadfitting.object=struct('Style','pushbutton','String','load');
pard.loadfitting.position=[8,1];
pard.loadfitting.Width=1;
pard.loadfitting.Height=1;
pard.loadfitting.tab='tab1';

pard.savefitting.object=struct('Style','pushbutton','String','save');
pard.savefitting.position=[8,2];
pard.savefitting.Width=1;
pard.savefitting.Height=1;
pard.savefitting.tab='tab1';

end

function pard=guidefmodel(obj,number)
pard.tab.tabmodel='Model';
pard.modelname.object=struct('Style','edit','String','');
pard.modelname.position=[4,1];
pard.modelname.Width=2;
pard.modelname.tab=['tabmodel_' num2str(number)];

pard.modelload.object=struct('Style','pushbutton','String','load model','Callback',{{@loadmodel_callback,obj}});
pard.modelload.position=[4,3];
pard.modelload.Width=1.3;
pard.modelload.tab=['tabmodel_' num2str(number)];

pard.layert.object=struct('Style','text','String','layer');
pard.layert.position=[5,1];
pard.layert.Width=2;
pard.layert.tab=['tabmodel_' num2str(number)];

pard.layer.object=struct('Style','edit','String','1');
pard.layer.position=[5,3];
pard.layer.Width=1;
pard.layer.tab=['tabmodel_' num2str(number)];

pard.pixelsizefitt.object=struct('Style','text','String','Pixel size');
pard.pixelsizefitt.position=[6,1];
pard.pixelsizefitt.Width=2;
pard.pixelsizefitt.tab=['tabmodel_' num2str(number)];

pard.pixelsizefit.object=struct('Style','edit','String','10');
pard.pixelsizefit.position=[6,3];
pard.pixelsizefit.Width=1;
pard.pixelsizefit.tab=['tabmodel_' num2str(number)];

pard.eps_fitt.object=struct('Style','text','String','eps');
pard.eps_fitt.position=[7,1];
pard.eps_fitt.Width=1;
pard.eps_fitt.tab=['tabmodel_' num2str(number)];

pard.eps_fit.object=struct('Style','edit','String','5e-6');
pard.eps_fit.position=[7,3];
pard.eps_fit.Width=1;
pard.eps_fit.tab=['tabmodel_' num2str(number)];

pard.weight_fitt.object=struct('Style','text','String','Weight');
pard.weight_fitt.position=[8,1];
pard.weight_fitt.Width=1;
pard.weight_fitt.tab=['tabmodel_' num2str(number)];

pard.weight_fit.object=struct('Style','edit','String','1');
pard.weight_fit.position=[8,3];
pard.weight_fit.Width=1;
pard.weight_fit.tab=['tabmodel_' num2str(number)];

pard.tab.tabpar='Parameters';
pard.tab.tabconv='Convert';
end