classdef WeiMoFit<interfaces.SEEvaluationProcessor
%     Fits circular model to localizations in the ROI to determine the
%     radius and center coordinates of circular structures.
    properties
        pe
        locmopath
        api
    end
    methods
        function obj=WeiMoFit(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});

                %set paths, initialize according to run.m
                obj.locmopath=rel2abs("../LocMoFitPy2/");
                addpath(obj.locmopath);
                try 
                    obj.pe = pyenv(Version="../LocMoFitPy2/.pixi/envs/default/bin/python", ExecutionMode= "OutOfProcess");
                catch err
                    terminate(pyenv)
                     obj.pe = pyenv(Version="../LocMoFitPy2/.pixi/envs/default/bin/python", ExecutionMode= "OutOfProcess");
                end
                if count(py.sys.path, obj.locmopath) == 0
                    insert(py.sys.path, int32(0), obj.locmopath);
                end
                obj.api = py.importlib.import_module("matlab_api");

        end
        function out=run(obj, p)
           out=runintern(obj,p);
           out=[];
        end
     
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end

end

function  out=runintern(obj,p)
layers=find(obj.getPar('sr_layerson'));
locs=obj.getLocs({"xnm","ynm","znm","xnmerr","ynmerr","locprecznm"},'layer',layers,'size',p.se_siteroi/2);  

mpos=horzcat(mean(locs.xnm),mean(locs.ynm),mean(locs.znm));

locarray=horzcat(locs.xnm-mpos(1),locs.ynm-mpos(2),locs.znm-mpos(3));
loc_precisions=horzcat(locs.xnmerr,locs.ynmerr,locs.locprecznm);
kwargs = py.dict;
% kwargs{"C"}=1e-6;

freeze = py.tuple({});
tic
res = obj.api.run_locmofit( ...
    "SphericalCap", ...
    locarray, ...
    loc_precisions, ...
    pyargs( ...
    "seed", int64(3), ...
    "model_init_kwargs", kwargs, ...
    "freeze", freeze, ...
    "max_iter", int64(200), ...
    "tol", 1e-12 ...
    ) ...
    );
toc

model=double(res{"positions"});
if obj.display
    ax= obj.setoutput('fit');
    hold(ax,'off')
    plot(ax,locs.xnm-mpos(1),locs.ynm-mpos(2),'k.');
    hold(ax,'on')
    ch=convhull(model(:,1:2));
    plot(ax,model(ch,1),model(ch,2))
    plot(ax,model(:,1),model(:,2),'m.')
    
    % fitpos=double(res{})
    
    plot(ax,locs.xnm-mpos(1),locs.znm-mpos(3)-p.se_siteroi,'k.');
    plot(ax,model(:,1),model(:,3)-p.se_siteroi,'m.')

    plot(ax,locs.ynm-mpos(2)-p.se_siteroi,locs.znm-mpos(3)-p.se_siteroi,'k.');
    plot(ax,model(:,2)-p.se_siteroi,model(:,3)-p.se_siteroi,'m.')

end
out=res;

end


function absPath = rel2abs(rel)
try
cur = pwd;
c = onCleanup(@() cd(cur));          % ensure we return to original folder
cd(rel);
absPath = pwd;
catch err
    absPath="";
end
end


function pard=guidef(obj)

pard.ch1.object=struct('Style','checkbox','String','checkbox example','Value',1);
pard.ch1.position=[2,1];
pard.ch1.Width=2;

pard.t1.object=struct('Style','text','String','text example');
pard.t1.position=[3,1];
pard.t1.Width=2;

pard.v1.object=struct('Style','edit','String','0');
pard.v1.position=[3,3];
pard.v1.Width=1;

pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
pard.plugininfo.type='ROI_Evaluate';
pard.plugininfo.description='';
end
