classdef CorrectDepthDependentOffset<interfaces.DialogProcessor&interfaces.SEProcessor
    % CorrectDepthDependentOffset corrects the z coordinates of
    % localizations based on Nup96's ring separation, 50 nm.
    % This plugin requires SMLMModelFitGUI been run.
    methods
        function obj=CorrectDepthDependentOffset(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
%             obj.showresults=false;
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','CorrectDepthDependentOffset');
            notify(obj.P,'backup4undo');
            
            %% Get info from SMLMModelFit
            sites = obj.SE.sites;
            lUsed = getFieldAsVector(sites, 'annotation.use');
            usedSites = sites(lUsed);
            
            evalList = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
            indProcessor = find(strcmp('SMLMModelFitGUI',evalList));
            fitter = obj.SE.processors.eval.processors{indProcessor}.fitter;
            
            [~,idxRingDist] = fitter.wherePar('pars.m2.lPar.z');
            ringDistS1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxRingDist);
            
            % get z pos of the 1st model
            [~,idxZ] = fitter.wherePar('pars.m1.lPar.z');
            zM1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxZ);
            
            % move the z pos to the center of the two rings
            z = zeros(size(zM1));
            for k = length(usedSites):-1:1
                fitter.allParsArg = usedSites(k).evaluation.SMLMModelFitGUI.allParsArg;
                z(k) = zM1(k)+fitter.rel(fitter.wherePar('pars.m2.lPar.z'),3,3)/2;
            end
            
            %% determine the cutoff of one-ring NPCs
            ax=obj.initaxis('One-ring NPC');
            [intersection, bin_edges,~,curve] = getOneRingCutoff(ringDistS1);
            histogram(ax,ringDistS1,bin_edges)
            hold(ax,'on')
            plot(curve)
            h = xline(intersection);
            legend(h,'cutoff')
            hold(ax,'off')
            xlabel('ring separation (nm)')
            ylabel('count')
            title(['Cutoff: ' sprintf('%2.1f', intersection)])
            
            lTwoRing = ringDistS1>intersection;

            %% moving mean
            xx = -120:10:80;
            xxCenter = movmean(xx,2);
            dz = ((ringDistS1(lTwoRing))-50)/50;
            yy=bindata(z(lTwoRing),dz,xx,'mean');
            
            ax2=obj.initaxis('Depth-depedent offset');
            
            plot(ax2,z(lTwoRing),dz, ' ob')
            hold(ax2, 'on')
            h1 = plot(ax2, xx,yy);
            hold(ax2, 'off')
            
            
            %% curve fit
            f = fit(z(lTwoRing),dz,'poly1','Robust','LAR');
            hold(ax2, 'on')
            h2 = plot(f, 'c');
            hold(ax2, 'off')
            xlabel(ax2, 'z position (nm)')
            ylabel(ax2, 'ring separation (nm)')
            legend([h1 h2],{'Move mean','Linear fit'})
            
            z0 = (0-f.p2)/f.p1;
            fInt = polyint([f.p1 f.p2]);
            c = -(fInt(1)*z0^2+fInt(2)*z0);
            fInt(3) = c;

            dzFun = @(z) fInt(1).*z.^2+fInt(2).*z+fInt(3);
            
            
            %% apply the correction
            if ~p.preview
                locs = obj.locData.loc;
                obj.locData.loc.znm = locs.znm-dzFun(locs.znm);           
                obj.locData.regroup;
                obj.locData.filter;
            end
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef

pard.preview.object=struct('String','Preview','Style','checkbox');
pard.preview.position=[2,1];
pard.preview.Width=1.5;

pard.plugininfo.type='ProcessorPlugin';


end