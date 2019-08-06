classdef LocIndices<interfaces.SEEvaluationProcessor
%     Reads out the indices of the localizations inside the ROI and stores
%     them. This can greatly accelarate evaluation of data sets with a very
%     large number of localizations.
    properties
        
    end
    methods
        function obj=LocIndices(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function makeGui(obj)
            makeGui@interfaces.SEEvaluationProcessor(obj);
%             obj.guihandles.saveimagesb.Callback={@saveimagesb_callback,obj};
        end
        function out=run(obj,p)
            try
            out=runintern(obj,p);
            catch err
                err
                out=[];
            end
         
        end
        
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef

pard.roiselect.object=struct('Style','popupmenu','String',{{'site roi round','site roi square','site fov','freeroi'}},'Value',2);
pard.roiselect.position=[1,1];
pard.roiselect.Width=2;

pard.plugininfo.type='ROI_Evaluate';
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
pard.plugininfo.description='Reads out the indices of the localizations inside the ROI and stores them. This can greatly accelarate evaluation of data sets with a very large number of localizations.';
end


function out=runintern(obj,p)
layers=find(obj.getPar('sr_layerson'));
switch p.roiselect.Value
    case 1 %square roi
        sizeh=[1 1]*p.se_siteroi;
    case 2 %round roi
        sizeh=p.se_siteroi;
    case 3 %FoV*
        sizeh=[1 1]*p.se_sitefov;
    case 4 %freeroi
        sizeh='freeroi';
        
end

[locs]=obj.getLocs({'ingrouped','inungrouped'},'layer',layers,'size',sizeh);
out.ingrouped=find(locs.ingrouped);
out.inungrouped=find(locs.inungrouped);

end

