classdef ParameterInterface<handle
    %global parameter handling
    properties
        P %interfaces.ParameterData. Parameters are saved here
        inputParameters %parameters a module needs from P. Passed on to getAllParameters etc
        outputParameters %parameter a module provides. 
    end
    methods
        function attachPar(obj,par)
            %attaches parameter object to module
            obj.P=par;
        end
        function setPar(obj,field,value)
            %sets global parameter
            %setPar(parameter,value)
            obj.P.par.(field).content=value;
            obj.P.par.(field).isGuiPar=false;
            obj.P.par.(field).synchronize=false;
            obj.P.par.(field).isGuiPar=0;
            obj.P.par.(field).synchronize=0;
            obj.P.par.(field).handle=[];
            obj.P.par.(field).syncmode=[];
            obj.P.par.(field).obj=obj;
            obj.P.par.(field).changecallback={};
        end
        function value=getPar(obj,field,tested)
            %gets global parameter
            %getPar(obj,field)           
            if nargin<3
                tested=false;
            end
            if tested||isfield(obj.P.par,field)
                value=obj.P.par.(field).content;
            else
                value=[];
            end
        end
        function p=getAllParameters(obj,inputParameters)   
            %p contains all parameters described by inputParaemters. parameter names as fields
            % p=getAllParameters(inputParameters). 
            %inputParameters (cell array) optional. If empty, use obj.inputParameters
            if nargin<2
                inputParameters=obj.inputParameters;
            end
            p=[];
            for k=1:length(inputParameters)
                p.(inputParameters{k})=obj.getPar(inputParameters{k});
            end
        end
        function updateObjectParameters(obj)
            %writes obj.outputParameters to Parameters
            %TODO rarely used. Take out?
            for k=1:length(obj.outputParameters)
                if isprop(obj,obj.outputParameters{k})
                    obj.setPar(obj.outputParameters{k},obj.(obj.outputParameters{k}));
                end
            end
        end
        
        function addResults(obj,field,value)
            if isfield(obj.P,field)
                obj.P.results.(field)(end+1)=value;
            else
                obj.setResults(field,value)
            end
        end
        
        function setResults(obj,field,value)
            obj.P.results.(field)=value;
        end
        function value=getResults(obj,field)
            if isfield(obj.P.results,field)
            value=obj.P.results.(field);
            else
                value=[];
            end
        end
        
        function setAutoResults(obj,pluginpath,value)
            q=value;
           
            for k=length(pluginpath):-1:2
                 p=[];
                p.(pluginpath{k})=q;
                q=p;
            end
            obj.P.autoresults.(pluginpath{1})=p;
        end
       
    end
end
        