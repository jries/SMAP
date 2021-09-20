classdef RemoveFields<interfaces.DialogProcessor
%     Define any fields as main x, y and z coordinates
    methods
        function obj=RemoveFields(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p) 
            out=[];
            obj.setPar('undoModule','Mathematics');
            notify(obj.P,'backup4undo');
            changes={'changex','changey','changez'};
            field2Rm = {};
            for k = 1:5
                if p.(['rmField' num2str(k)])
                	field2Rm = [field2Rm p.(['field' num2str(k)]).selection];
                end
            end
            obj.locData.loc = rmfield(obj.locData.loc, field2Rm);
            obj.locData.grouploc = rmfield(obj.locData.grouploc, field2Rm);
            obj.setPar('locFields',fieldnames(obj.locData.loc));
            obj.locData.regroup;
            obj.locData.filter;
                
        end
        function initGuiFinal(obj)  
%             fields={'rmField1','rmField2','rmField3','rmField4','rmField5'};
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard=guidef(obj)

pard.rmField1.object=struct('String','','Style','checkbox');
pard.rmField1.position=[1,1];

pard.field1.object=struct('Style','popupmenu','String',' ');
pard.field1.position=[1,2];

pard.rmField2.object=struct('String','','Style','checkbox');
pard.rmField2.position=[2,1];

pard.field2.object=struct('Style','popupmenu','String',' ');
pard.field2.position=[2,2];

pard.rmField3.object=struct('String','','Style','checkbox');
pard.rmField3.position=[3,1];

pard.field3.object=struct('Style','popupmenu','String',' ');
pard.field3.position=[3,2];

pard.rmField4.object=struct('String','','Style','checkbox');
pard.rmField4.position=[4,1];

pard.field4.object=struct('Style','popupmenu','String',' ');
pard.field4.position=[4,2];

pard.rmField5.object=struct('String','','Style','checkbox');
pard.rmField5.position=[5,1];

pard.field5.object=struct('Style','popupmenu','String',' ');
pard.field5.position=[5,2];

pard.syncParameters={{'locFields','field1',{'String'}},...
    {'locFields','field2',{'String'}},...
    {'locFields','field3',{'String'}},...
    {'locFields','field4',{'String'}},...
    {'locFields','field5',{'String'}}};


pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Define any fields as main x, y and z coordinates';
end