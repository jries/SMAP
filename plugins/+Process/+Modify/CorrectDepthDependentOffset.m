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
            
            knownSeparation = 49.3;
            
            %% Get info from SMLMModelFit
            se = obj.SE;
            sites = obj.SE.sites;
            lUsed = getFieldAsVector(sites, 'annotation.use');
            usedSites = sites(lUsed);
            
            evalList = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
            indProcessor = find(strcmp('SMLMModelFitGUI',evalList));
            fitter = obj.SE.processors.eval.processors{indProcessor}.fitter;
            
            ID = getFieldAsVector(usedSites, 'ID');
            
            [~,idxRingDist] = fitter.wherePar('pars.m2.lPar.z');
            ringDistS1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxRingDist);
            
            [~,idxXrot] = fitter.wherePar('pars.m1.lPar.xrot');
            xrot = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxXrot);
            
            [~,idxYrot] = fitter.wherePar('pars.m1.lPar.yrot');
            yrot = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxYrot);
            
            % find the true elevation
            [x,y,z] = rotcoord3(0,0,-1, deg2rad(xrot), deg2rad(yrot), 0, 'XYZ');
            [~,eleOri,~] = cart2sph(0,0,-1);
            [~,ele,~] = cart2sph(x,y,z);
            trueEle = -(ele-eleOri);
            
            ringDistS1Z = ringDistS1.*cos(trueEle');
            knownSepZ = knownSeparation.*cos(trueEle');
            
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
            xx = getHistogramEdge(z(lTwoRing),10);
            xxCenter = movmean(xx,2);
            dOffset = (ringDistS1Z-knownSepZ)./ringDistS1Z;
            dOffset = dOffset(lTwoRing);
            yy=bindata(z(lTwoRing),dOffset,xx,'mean');
            
            ax2=obj.initaxis('Depth-depedent offset');
            
            z_twoRing = z(lTwoRing);
            ID_twiRing = ID(lTwoRing);
            lFocus = z_twoRing<100&z_twoRing>-100;
%             plot(ax2,z_twoRing(lFocus),dOffset(lFocus), ' ob')
            plotSElink(ax2,z_twoRing(lFocus),dOffset(lFocus),ID_twiRing(lFocus),se, ' ob');
            hold(ax2, 'on')
            plot(ax2,z_twoRing(~lFocus),dOffset(~lFocus), ' ok')
            h1 = plot(ax2, xx,yy);
            hold(ax2, 'off')
            
            
            %% curve fit
            f = fit(z_twoRing(lFocus),dOffset(lFocus),'poly1','Robust','LAR');
            hold(ax2, 'on')
            h2 = plot(f, 'c');
            hold(ax2, 'off')
            xlabel(ax2, 'z position (nm)')
            ylabel(ax2, 'Relative offset')
            legend([h1 h2],{'Move mean','Linear fit'})
            
            z0 = (0-f.p2)/f.p1;
            fInt = polyint([f.p1 f.p2]);
            c = -(fInt(1)*z0^2+fInt(2)*z0);
            fInt(3) = c;

            dOffsetFun = @(z) fInt(1).*z.^2+fInt(2).*z+fInt(3);
            
            %%
            ax2_1=obj.initaxis('Ring separation');
            plotSElink(ax2_1, z_twoRing,ringDistS1Z(lTwoRing),ID,se,' ob') 
            f = fit(z_twoRing,ringDistS1Z(lTwoRing),'poly1','Robust','LAR');
            hold(ax2_1, 'on')
            h3 = plot(f, 'c');
            hold(ax2_1, 'off')
            
            %% validation
            newZLower = zM1-dOffsetFun(zM1);
            
            zUpper = zM1+ringDistS1Z;
            newZUpper = zUpper-dOffsetFun(zUpper);
            newRingDistS1Z = newZUpper-newZLower;
            
            lateralDist = sqrt(ringDistS1.^2-ringDistS1Z.^2);
            newRingDist = sqrt(lateralDist.^2+newRingDistS1Z.^2);
            newZPos = newZLower+newRingDistS1Z/2;
            f2 = fit(newZPos(lTwoRing),newRingDist(lTwoRing),'poly1','Robust','LAR');
            ax3=obj.initaxis('Validation');
            plot(ax3, newZPos(lTwoRing),newRingDist(lTwoRing),' ob')
            hold(ax3, 'on')
            h2 = plot(f2, 'c');
            hold(ax3, 'off')
            xlabel(ax3, 'z position (nm)')
            ylabel(ax3, 'ring separation (nm)')
            legend([h1 h2],{'Move mean','Linear fit'})
            
            
            
            %% apply the correction
            if ~p.preview
                obj.setPar('undoModule','CorrectDepthDependentOffset');
                notify(obj.P,'backup4undo');
            
                locs = obj.locData.loc;
                obj.locData.loc.znm = locs.znm-dOffsetFun(locs.znm);           
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