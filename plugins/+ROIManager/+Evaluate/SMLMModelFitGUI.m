classdef SMLMModelFitGUI<interfaces.SEEvaluationProcessor

    methods
        function obj=SMLMModelFitGUI(varargin)
            obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        
        function out=run(obj, inp, varargin)
            p = inputParser;
            addParameter(p,'onlySetUp',false, @islocgical);
            parse(p,varargin{:});
            results = p.Results;
            out=runSMLMModelFitGUI(obj, inp, results.onlySetUp);
        end
        
        function makeGui(obj,varargin)
            obj.setPar('currentLoadedModel',[]);
            if nargin >2 && varargin{2}==true %gui of tabs
                 makeGui@interfaces.GuiModuleInterface(obj,varargin{1});
            else %real make GUI of main gui
                makeGui@interfaces.GuiModuleInterface(obj); %make the main GUI from the guidef definitions
                %Settings               
                optimizernames={'fminsearch','particleswarm'};%if possible, get from fitter object
                obj.guihandles.optimizer.String=optimizernames;
                Data={'Display','iter';'UseVectorized','true'};
                oldh=obj.guihandles.optimizerpar;
                pos=oldh.Position;
                htable=uitable(oldh.Parent,'Data',Data,'Position',pos);
                htable.ColumnEditable = true;
                delete(oldh);
                obj.guihandles.optimizerpar=htable;
                
                % layer settings
                oldh = obj.guihandles.layerSetting;
                pos=oldh.Position;
                colNames={'Layer','Type','Label','Name','Fix','LB','UB','Min','Max','Value'};
                hLayer=uitable(oldh.Parent,'Data',{},'Position',pos);
                hLayer.ColumnName = colNames;
                hLayer.ColumnEditable = true;
                hLayer.CellEditCallback = {@layerSetting_callback,obj};
                delete(oldh);
                obj.guihandles.layerSetting=hLayer;

                %add first module and  + tab
                obj.guihandles.tabgroup=obj.guihandles.tab1.Parent;
                obj.addguitotab(1);
                obj.guihandles.addtab=uitab(obj.guihandles.tabgroup,'Title','+');
                obj.setPar('numMod',1);                 % init of the model counts
                obj.guihandles.tabgroup.SelectionChangedFcn={@selectLayer_callback,obj};
                
                %Convert               
                addconverttotab(obj);
                oldh=obj.guihandles.anchorConvert;
                pos=oldh.Position;
                htable=uitable(oldh.Parent,'Data',{},'Position',[pos(1:2)+[5 10] 300 200]);
                colNames={'Source', 'Rule', 'Target_fit', 'Target_usr'};
                htable.ColumnName = colNames;
                htable.CellEditCallback = {@convertTable_callback,obj};
                htable.ColumnFormat = {{'own'},[],{'none'},[]};
                htable.ColumnEditable = true;
                delete(oldh);
                obj.guihandles.anchorConvert=htable;
                
                % create the SMLMModelFit
                % check the dim of the data and update the corresponding
                % seeting in the fitter
                if isfield(obj.locData.loc,'znm')
                    dataDim = 3;
                else
                    dataDim = 2;
                end
                fitter = SMLMModelFit('DataDim',dataDim);
                obj.setPar(['fitter_' obj.name], fitter);
            end
        end
        function pard=guidef(obj)
            pard=SMLMModelFitGUIdef(obj);
        end
        
        function addguitotab(obj,number)
            run_addguitotab(obj,number);
        end
        
        function addconverttotab(obj)
            run_addconverttotab(obj);
        end
    end

end

function selectLayer_callback(tabgroup,eventdata,obj) 
    % if + tab selected this makes a new model
    layertitle=(eventdata.NewValue.Title);
    if strcmp(layertitle,'+')
        numtabs=length(tabgroup.Children);
        obj.addguitotab(numtabs-2)
        obj.setPar('numMod',numtabs-2)   % save the number of model
        s=1:length(tabgroup.Children);
        % shift the order of table
        s(end-2)=s(end);
        s(end)=s(end)-2;
        s(end-1)=s(end);
        s(end)=s(end)+1;
        tabgroup.Children=tabgroup.Children(s);
        tabgroup.SelectedTab=tabgroup.Children(end-2); 
    end
end

function convertTable_callback(a,b,obj) 
    if b.Indices(2) == 4
        htable = obj.guihandles.anchorConvert;
        optionTarget = unique([htable.ColumnFormat{3} {['usr_' b.NewData]}]);
        htable.ColumnFormat{3} = optionTarget;
        obj.guihandles.anchorConvert=htable;
    end
end

function layerSetting_callback(a,b,obj) 
    fitter = obj.getPar(['fitter_' obj.name]);
    fn={'model','type','label','name','fix','lb','ub','min','max','value'};
    indices = b.Indices;
    layer = a.Data{indices(1),1};
    layer = strrep(layer,'layer','9');
    [~,idx] = fitter.wherePar(['pars.m' layer '.offset.weight']);
    fitter.allParsArg.(fn{indices(2)})(idx) = b.NewData;
    obj.setPar(['fitter_' obj.name], fitter);
end