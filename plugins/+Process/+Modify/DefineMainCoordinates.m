classdef DefineMainCoordinates<interfaces.DialogProcessor
%     Define any fields as main x, y and z coordinates
    methods
        function obj=DefineMainCoordinates(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p) 
            out=[];
            obj.setPar('undoModule','Mathematics');
            notify(obj.P,'backup4undo');
            changes={'changex','changey','changez'};
            targetfields={'xnm','ynm','znm'};
            targetprec={'locprecnm','locprecnm','locprecznm'};
            fields={'xfield','yfield','zfield'};
            errfields={'xerrfield','xerrfield','zerrfield'};
            facfields={'xfac','yfac','zfac'};
            
            newfields={p.xfield.selection,p.yfield.selection,p.zfield.selection,p.xerrfield.selection,p.zerrfield.selection};
            
            loctemp=copyfields([],obj.locData.loc,targetfields);
            loctemp=copyfields(loctemp,obj.locData.loc,targetprec);
             loctemp=copyfields(loctemp,obj.locData.loc,newfields);           
            for k=1:3
                if p.(changes{k})
                    if p.changefield
                        backupfield=[targetfields{k} '_original'];
                        if ~isfield(obj.locData.loc,backupfield) && isfield(obj.locData.loc,targetfields{k})
                            obj.locData.setloc(backupfield,loctemp.(targetfields{k}));
                        end

                        fac=p.(facfields{k});
                        if length(fac)>1
                            off=fac(2);
                            fac=fac(1);
                        else
                            off=0;
                        end
                        xnew=(loctemp.(p.(fields{k}).selection)-off)*fac;
                        obj.locData.setloc(targetfields{k},xnew);
                    end
                    if p.changeprec
                        backupfield=[targetprec{k} '_original'];
                        if ~isfield(obj.locData.loc,backupfield) && isfield(obj.locData.loc,targetfields{k})
                            obj.locData.setloc(backupfield,loctemp.(targetprec{k}));
                        end
                        fac=p.(facfields{k})(1);
                        xnew=(loctemp.(p.(errfields{k}).selection))*fac;
                        obj.locData.setloc(targetprec{k},xnew);
                    end
                end
            end
            obj.locData.regroup;
            obj.locData.filter;
        end
        function initGuiFinal(obj)  
            fields={'xfield','yfield','zfield','xerrfield','zerrfield'};
            fieldnames={'xnm','ynm','znm','locprecnm','locprecznm'};
            for k=1:length(fields)
                stri1=obj.guihandles.(fields{k}).String;
                indg=find(contains(stri1,fieldnames{k}));
                if ~isempty(indg)
                    obj.guihandles.(fields{k}).Value=indg(1);
                end
            end   
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function revert_callback(obj,a,b)
            targetfields={'xnm','ynm','znm','locprecnm','locrpecznm'};
            for k=1:length(targetfields)
                backupfield=[targetfields{k} '_original'];
                if isfield(obj.locData.loc,backupfield)
                    obj.locData.setloc(targetfields{k},obj.locData.loc.(backupfield));
                end
            end
            obj.locData.regroup;
            obj.locData.filter;
        end
    end
end

function pard=guidef(obj)

pard.changex.object=struct('String','xnm ','Style','checkbox');
pard.changex.position=[2,2];
pard.changey.object=struct('String','ynm ','Style','checkbox');
pard.changey.position=[2,3];
pard.changez.object=struct('String','znm ','Style','checkbox');
pard.changez.position=[2,4];

pard.changefield.object=struct('String','field','Style','checkbox');
pard.changefield.position=[3,1];

pard.xfield.object=struct('Style','popupmenu','String',' ');
pard.xfield.position=[3,2];
pard.yfield.object=struct('Style','popupmenu','String',' ');
pard.yfield.position=[3,3];
pard.zfield.object=struct('Style','popupmenu','String',' ');
pard.zfield.position=[3,4];

pard.changeprec.object=struct('String','precision','Style','checkbox');
pard.changeprec.position=[4,1];
pard.xerrfield.object=struct('Style','popupmenu','String',' ');
pard.xerrfield.position=[4,2];
pard.xerrfield.Width=2;
pard.zerrfield.object=struct('Style','popupmenu','String',' ');
pard.zerrfield.position=[4,4];


pard.factext.object=struct('String','factor (,offset) ','Style','text');
pard.factext.position=[5,1];

pard.xfac.object=struct('String','1 0','Style','edit');
pard.xfac.position=[5,2];

pard.yfac.object=struct('String','1 0','Style','edit');
pard.yfac.position=[5,3];
pard.zfac.object=struct('String','1 0','Style','edit');
pard.zfac.position=[5,4];


pard.revert.object=struct('String','Revert to original','Style','pushbutton','Callback',@obj.revert_callback);
pard.revert.position=[6,2];
pard.revert.Width=2;
pard.syncParameters={{'locFields','xfield',{'String'}},...
    {'locFields','yfield',{'String'}},...
    {'locFields','zfield',{'String'}},...
    {'locFields','xerrfield',{'String'}},...
    {'locFields','zerrfield',{'String'}}};


pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Define any fields as main x, y and z coordinates';
end