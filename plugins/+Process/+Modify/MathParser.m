classdef MathParser<interfaces.DialogProcessor
    properties
        equationhistory
        historyfile='settings/temp/MathParser.txt';
    end
    methods
        function obj=MathParser(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        function initGui(obj)
            if exist(obj.historyfile,'file')
                obj.equationhistory=readtable(obj.historyfile);
                p=obj.getGuiParameters;
                p.resultfieldh.String=obj.equationhistory.resultfield;
                p.equationh.String=obj.equationhistory.equation;
                obj.setGuiParameters(p);
            else
                tt=struct('resultfield',{{''}},'equation',{{''}});
                obj.equationhistory=struct2table(tt);
            end
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','Math Parser');
            notify(obj.P,'backup4undo');
            locs = obj.locData.loc;
            evalstr=p.equation;
            fn=fieldnames(locs);
            
            if p.dataselect_all
                indin=true(size(locs.filenumber));
            else
                indin=locs.filenumber==p.dataselect.Value;
            end
            nl=cellfun(@length,fn);
            [~ ,sortind]=sort(nl,'descend');
            fn=fn(sortind);
            for k = 1:length (fn)
                ind=strfind(evalstr,fn{k});
                for m=length(ind):-1:1
                    if ind(m)>1
                        lbefore=evalstr(ind(m)-1);
                    else
                        lbefore='';
                    end
                    indend=ind(m)+length(fn{k});
                    if indend<length(evalstr)
                        lafter=evalstr(indend+1);
                    else
                        lafter='';
                    end
                    
                    if strcmp(lbefore,'?') || strcmp(lafter,'?') 
                        continue
                    end
                    if isempty(lbefore) || (lbefore~='.' && ~isstrprop(lbefore,'alpha')) 
                    evalstr=[evalstr(1:ind(m)-1) 'locs.?' evalstr(ind(m):ind(m)+length(fn{k})-1) '?(indin)' evalstr(ind(m)+length(fn{k}):end)];
                    end
%                     evalstr = strrep(evalstr,[lbefore fn{k}],[lbefore 'locs.'  fn{k}]);
                end
            end
            evalstr(evalstr=='?')=[];
             try
                newval=eval(evalstr);
                if isfield(obj.locData.loc,p.resultfield)
                    obj.locData.loc.(p.resultfield)(indin)=newval;
                else
                 obj.locData.setloc(p.resultfield,newval,indin);
                end
                 obj.locData.filter(p.resultfield)
                 obj.locData.regroup;
                 obj.setPar('locFields',fieldnames(obj.locData.loc))
               
                 exe=(contains(obj.equationhistory.equation,p.equation));
%                  exr=find(obj.equationhistory.resultfield,p.resultfield);
                 if sum(exe)==0
                     l=length(obj.equationhistory.equation);
                     obj.equationhistory.resultfield(2:min(l,9)+1)=obj.equationhistory.resultfield(1:min(l,9));
                     obj.equationhistory.equation(2:min(l,9)+1)=obj.equationhistory.equation(1:min(l,9));
                     obj.equationhistory.resultfield{1}=p.resultfield;
                     obj.equationhistory.equation{1}=p.equation;
                 else
                     
                     obj.equationhistory.resultfield(2:sum(~exe)+1)=obj.equationhistory.resultfield(~exe);
                     obj.equationhistory.equation(2:sum(~exe)+1)=obj.equationhistory.equation(~exe);
                     obj.equationhistory.resultfield{1}=p.resultfield;
                     obj.equationhistory.equation{1}=p.equation;
                     
%                      obj.equationhistory.resultfield(exe)=p.resultfield;
                 end
                 
                 writetable(obj.equationhistory,obj.historyfile);
                 p.resultfieldh.String=obj.equationhistory.resultfield;
                 p.equationh.String=obj.equationhistory.equation;
                 p.resultfieldh.Value=1;
                 p.equationh.Value=1;
                 obj.setGuiParameters(p);
                 
             catch err
                 disp(['could not evaluate equation: ' evalstr])
                 disp(err.message)
             end
             
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function dobutton(obj,callobj,event,inp1)
            disp('dobutton')
           x= obj.locData.loc.xnm;
           mean(x)
        end
        function history_callback(obj,a,b)
            neval=a.Value;
            p=obj.getGuiParameters;
            p.resultfield=p.resultfieldh.String{neval};
            p.equation=p.equationh.String{neval};
            obj.setGuiParameters(p);
        end
    end
end


        

function pard=guidef(obj)
pard.t1.object=struct('String','result field','Style','text');
pard.t1.position=[2,1];

pard.resultfield.object=struct('String','within','Style','edit');
pard.resultfield.position=[3,1];
pard.resultfield.Width=0.8;


pard.t2.object = struct('String','=','Style','text');
pard.t2.position=[3,1.8];
pard.t2.Width=0.2;

pard.t3.object=struct('String','Equation. Use fieldnames (e.g. xnm, phot) as variables','Style','text');
pard.t3.position=[2,2];
pard.t3.Width=3;

pard.equation.object = struct('String','(locprecnm<25 & PSFxnm>100) | numberInGroup>1','Style','edit');
pard.equation.position=[3,2];
pard.equation.Width=3;


pard.resultfieldh.object=struct('String',{{''}},'Style','popupmenu','Callback',@obj.history_callback);
pard.resultfieldh.position=[4,1];
pard.resultfieldh.Width=0.8;


pard.equationh.object = struct('String',{{''}},'Style','popupmenu','Callback',@obj.history_callback);
pard.equationh.position=[4,2];
pard.equationh.Width=3;


pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[1,1];
pard.dataselect.object.TooltipString='choose localization file data set';

pard.dataselect_all.object=struct('Style','checkbox','String','all');
pard.dataselect_all.position=[1,2];
pard.dataselect_all.object.TooltipString='choose localization file data set';

pard.syncParameters={{'filelist_short','dataselect',{'String'}},{'MathParserHistory','resultfieldh',{'Value'}},{'MathParserHistory','equationh',{'Value'}}};

pard.plugininfo.type='ProcessorPlugin';
end