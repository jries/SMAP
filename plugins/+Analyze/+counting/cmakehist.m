classdef cmakehist<interfaces.DialogProcessor
    methods
        function obj=cmakehist(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
                obj.inputParameters={'mainfile'};
        end
        
        function out=run(obj,p)
            clusters=obj.getResults('counting_clusters');
            histogram=cluster_mmaple_makehist(p,clusters);
            obj.setResults('counting_histogram',histogram);
            out=clusters;

        end
        
        function refit_callback(obj)
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        
        function savehist_callback(obj,a,b)
            p=obj.getAllParameters;
            histogram=obj.getResults('counting_histogram');
            of=p.mainfile;
            if isempty(of)
                off='*_hist.mat';
            else
                [path,f,ext]=fileparts(of);
                off=[path filesep f '_hist.mat'];
            end
            [f,path]=uiputfile(off);
            if f
                save([path,f],'histogram')
            end
            
        end
    end
end




function pard=guidef(obj)

pard.text1.object=struct('String','PSFsigma (nm) min:max','Style','text');
pard.text1.position=[1,1];
pard.text1.Width=2;
pard.c_PSFmin.object=struct('String','100','Style','edit');
pard.c_PSFmin.position=[1,3];
pard.c_PSFmin.Width=0.5;
pard.c_PSFmax.object=struct('String','130','Style','edit');
pard.c_PSFmax.position=[1,3.5];
pard.c_PSFmax.Width=0.5;


pard.text2.object=struct('String','Std cluster size (nm) min:max','Style','text');
pard.text2.position=[2,1];
pard.text2.Width=2;
pard.c_stdymin.object=struct('String','0','Style','edit');
pard.c_stdymin.position=[2,3];
pard.c_stdymin.Width=0.5;
pard.c_stdymax.object=struct('String','250','Style','edit');
pard.c_stdymax.position=[2,3.5];
pard.c_stdymax.Width=0.5;

% pard.N0_v.object=struct('String','10','Style','edit');
% pard.N0_v.position=[2,3];
pard.c_groupfield.object=struct('String',{{'ungrouped','grouped','grouped max dist','blink remove'}},'Style','popupmenu');
pard.c_groupfield.position=[3,1];
pard.c_groupfield.Width=2;

pard.text3.object=struct('String','max dark time (frames)','Style','text');
pard.text3.position=[4,1];
pard.text3.Width=2;

pard.dftime.object=struct('String','20','Style','edit');
pard.dftime.position=[4,3];

pard.savehist.object=struct('String','Save histogram','Style','pushbutton','Callback',{{@obj.savehist_callback}});
pard.savehist.position=[5,1];

pard.plugininfo.name='make brightness histogram';
pard.plugininfo.type='ProcessorPlugin';
end