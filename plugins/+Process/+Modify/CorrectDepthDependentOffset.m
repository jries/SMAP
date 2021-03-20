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
            %% internal parameters
            expectSep = 49.3;
            fixTwoRingCutoff = true;    % true for fixing the cutoff at twoRingCutoff; false for calculating the cutoff using getOneRingCutoff().
            twoRingCutoff = 35;
            
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
            
            % convert the rotation to the same as in the spherical
            % coordinate system
            [x,y,z] = rotcoord3(0,0,-1, deg2rad(xrot), deg2rad(yrot), 0, 'XYZ');
            [~,eleOri,~] = cart2sph(0,0,-1);
            [azi,ele,~] = cart2sph(x,y,z);
            sphEle = -(ele-eleOri);
            
            % z component of the ring distance
            ringDistS1Z = ringDistS1.*cos(sphEle');
            expectSepZ = expectSep.*cos(sphEle');
            
            % get z pos of model 1
            [~,idxZ] = fitter.wherePar('pars.m1.lPar.z');
            zM1 = getFieldAsVectorInd(usedSites, 'evaluation.SMLMModelFitGUI.allParsArg.value',idxZ);
            zUpper = zM1+ringDistS1Z;
            lateralDist = sqrt(ringDistS1.^2-ringDistS1Z.^2);
            
            % move the z pos to the mutual center of both rings
            z = zeros(size(zM1));
            for k = length(usedSites):-1:1
                fitter.allParsArg = usedSites(k).evaluation.SMLMModelFitGUI.allParsArg;
                z(k) = zM1(k)+fitter.rel(fitter.wherePar('pars.m2.lPar.z'),3,3)/2;
            end
            
            %% determine the cutoff of one-ring NPCs
            [intersection, bin_edges,~,curve] = getOneRingCutoff(ringDistS1);
            if ~fixTwoRingCutoff
                twoRingCutoff = intersection;
            end
            lTwoRing = ringDistS1>twoRingCutoff;
            lGood = lTwoRing&abs(z)<150;
            
            %% measured relative offset
            relMeasuredOffset = (ringDistS1Z-expectSepZ)./ringDistS1Z;
            relMeasuredOffset_good = relMeasuredOffset(lGood);
            z_good = z(lGood);
            ID_good = ID(lGood);

            %% Ring separation
            axRS=obj.initaxis('Ring separation vs z');
            plotSElink(axRS, z,ringDistS1,ID,se,' ob')
            hold(axRS,'on')
            h = yline(axRS, twoRingCutoff);
            hold(axRS,'off')
            xlabel(axRS, 'z position (nm)')
            ylabel(axRS, 'Ring separation (nm)')
            legend(h,{'Cutoff'})
            
            %% Reletive z offset
            axRO=obj.initaxis('Reletive offset');
            
            xx = getHistogramEdge(z,10);
            yy=bindata(z,relMeasuredOffset*100,xx,'mean');
            yy2=bindata(z,relMeasuredOffset*100,xx,'median');
            yy3=bindata(z,relMeasuredOffset*100,xx,'robustmean');
            
            plot(axRO, z,relMeasuredOffset*100,' ob') 
%             plotSElink(axRO, z,relMeasuredOffset*100,ID,se,' ob') 
            axRO.YLim = [-100 100];
            hold(axRO,'on')
            h1 = plot(axRO, xx,yy);
            h2 = plot(axRO, xx,yy2);
            h3 = plot(axRO, xx,yy3);
            h = yline(axRO, 100*(twoRingCutoff-expectSep)/twoRingCutoff);
            hold(axRO,'off')
            xlabel(axRO, 'z position (nm)')
            ylabel(axRO, 'Relative offset (%)')
            legend([h h1 h3 h2],{'Cutoff','Move mean','Move robustMean','Move median'})
            
            %% Correction factor
            axCF=obj.initaxis('Correction factor');
            correctionFactor = expectSepZ./ringDistS1Z;
            
            f2_CF = fit(z_good, correctionFactor(lGood),'poly2','Robust','LAR');
            correctedZ_CF = polyint([f2_CF.p1 f2_CF.p2 f2_CF.p3]);
            p4 = -polyval(correctedZ_CF,0);
            correctedZ_CF(4) = p4;
            correctedZFun = @(z) polyval(correctedZ_CF,z);
            
            pos_Lower = correctedZFun(zM1);
            pos_Upper = correctedZFun(zUpper);
            newRingDistS1Z_J = pos_Upper-pos_Lower;
            newRingDist_J = sqrt(lateralDist.^2+newRingDistS1Z_J.^2);
            
            yy=bindata(z(lGood),correctionFactor(lGood),xx,'median');
            
            h_TwoRing = plot(axCF, z(lGood),correctionFactor(lGood),' ob');
%             plotSElink(axRO, z,relMeasuredOffset*100,ID,se,' ob') 
            axCF.YLim = [0.5 3.5];
            hold(axCF,'on')
            h_NotTwoRing = plot(axCF, z(~lGood),correctionFactor(~lGood),' ok');
            h1 = plot(axCF, xx,yy);
            pf = plot(f2_CF,'- m');
            h = yline(axCF, expectSep/twoRingCutoff);
            hold(axCF,'off')
            xlabel(axCF, 'z position (nm)')
            ylabel(axCF, 'Relative offset (%)')
            legend([h_TwoRing h_NotTwoRing h h1 pf],{'Included','Excluded','Cutoff','Move median','Quadratic fit'})
            
            %% histogram
%             ax=obj.initaxis('Ring separation');
%             histogram(ax,ringDistS1,bin_edges)
%             hold(ax,'on')
%             plot(curve)
%             h = xline(twoRingCutoff);
%             legend(h,'cutoff')
%             hold(ax,'off')
%             xlabel('ring separation (nm)')
%             ylabel('count')
%             title(['Cutoff: ' sprintf('%2.1f', twoRingCutoff)])
%             
            %% Curve fit
            % moving mean
            xx = getHistogramEdge(z(lGood),10);
            yy=bindata(z(lGood),relMeasuredOffset_good,xx,'median');
            
            % Base plot
            %             lFocus = z_good<100&z_good>-100;
            ax2=obj.initaxis('Curve fit');
            h_good = plot(ax2,z(lGood),relMeasuredOffset(lGood), ' ob');
%             plotSElink(ax2,z_good(lFocus),relMeasuredOffset(lFocus),ID_twiRing(lFocus),se, ' ob');
            hold(ax2, 'on')
            h_NotTwoRing = plot(ax2, z(~lGood),relMeasuredOffset(~lGood),' ok');
%             plot(ax2,z_good(~lFocus),relMeasuredOffset(~lFocus), ' ok')
            h1 = plot(ax2, xx,yy);
            hold(ax2, 'off')
            
            % curve fit
            f = fit(z_good,relMeasuredOffset_good,'poly1','Robust','LAR');
            fP2 = fit(z_good,relMeasuredOffset_good,'poly2','Robust','LAR');
            hold(ax2, 'on')
            h3 = plot(fP2, '-m');
            hold(ax2, 'off')
            ax2.YLim = [-1 1];
            xlabel(ax2, 'z position (nm)')
            ylabel(ax2, 'Relative offset')
            legend([h_good h_NotTwoRing h1 h3],{'Included','Excluded','Move median', 'Quadratic fit'})
            
%             % first order
%             z0 = fzero(f, 0);
%             fInt = polyint([f.p1 f.p2]);
%             c = -(fInt(1)*z0^2+fInt(2)*z0);
%             fInt(3) = c;
%             relMeasuredOffsetFun = @(z) fInt(1).*z.^2+fInt(2).*z+fInt(3);
%             

            % 2nd order
            z0 = fzero(fP2, 0);
            if isnan(z0)
                z0 = fminsearch(@(x)abs(feval(fP2,x)),0);
            end
            fInt = polyint([fP2.p1 fP2.p2 fP2.p3]);
            c = -(fInt(1)*z0^3+fInt(2)*z0^2+fInt(3)*z0);
            fInt(4) = c;
            relMeasuredOffsetFun = @(z) polyval(fInt,z);
            
            %% Expected separation
            newZLower = zM1-relMeasuredOffsetFun(zM1);
            
            
            newZUpper = zUpper-relMeasuredOffsetFun(zUpper);
            newRingDistS1Z = newZUpper-newZLower;
            newRingDist = sqrt(lateralDist.^2+newRingDistS1Z.^2);
            newZPos = newZLower+newRingDistS1Z/2;
            f2 = fit(newZPos(lGood),newRingDist(lGood),'poly1','Robust','LAR');
            ax3=obj.initaxis('Expectation');
            p1 = plot(ax3, newZPos(lGood),newRingDist(lGood),' ob');
            hold(ax3, 'on')
            p2 = plot(ax3, newZPos(~lGood),newRingDist(~lGood),' ok');
            h2 = plot(f2, 'c');
            hold(ax3, 'off')
            xlabel(ax3, 'z position (nm)')
            ylabel(ax3, 'ring separation (nm)')
            legend([p1 p2 h2],{'Included','Discarded','Linear fit'})
            
            %% Jonas'
%             f_J = fit(z_good,expectSepZ(lGood)./ringDistS1Z(lGood),'poly1','Robust','LAR');
% %             hold(ax2_1, 'on')
% % %             h3 = plot(f_J, 'c');
% %             hold(ax2_1, 'off')
%             z0 = (0-f_J.p2)/f_J.p1;
%             fInt_J = polyint([f_J.p1 f_J.p2]);
%             c = 0;
%             fInt_J(3) = c;
           
            
            %% Comparison
            axComp=obj.initaxis('Method comparison');
            
            plot(axComp, newRingDist(lGood), newRingDist_J(lGood), ' ob')
            hold(axComp, 'on')
            leftCor = max(axComp.XLim(1),axComp.YLim(1));
            rightCor =  max(axComp.XLim(2),axComp.YLim(2));
            plot(axComp, [leftCor rightCor], [leftCor rightCor],'-k')
            hold(axComp, 'off')
            xlabel(axComp, 'Acts on offset')
            ylabel(axComp, 'Acts on znm')
            
%             %% derive model parameters
% 
%             newTrueEle = acos(newRingDistS1Z./newRingDist);
%             newEle = -newTrueEle+eleOri;
%             [v(:,1),v(:,2),v(:,3)] = sph2cart(azi',newEle,1);
% 
%             v_good = v(lGood,:);
%             rotAng = zeros([length(v_good) 2]);
%             xrot_good = xrot(lGood);
%             yrot_good = yrot(lGood);
%             axRot = obj.initaxis('rotation');
%             plot3(axRot, xrot(lGood),yrot(lGood), relMeasuredOffset_good,' ob')
%             for k = 1:length(v_good)
%                rotAng(k,:) = fminsearchbnd(@(x) rotConversion(x,v_good(k,:)),[xrot_good(k) yrot_good(k)],[-45 -45],[45 45]);
%             end
%             
%             sites_good = sites(lGood);
%             
%             newRingDist_good = newRingDist(lGood);
%             newZLower_good = newZLower(lGood);
            
            %% show the difference between measured and expectation
            if isfield(obj.locData.loc,'znm_rs_original')
                for k = length(usedSites):-1:1
                    expCorrRingSep(k) = usedSites(k).evaluation.CorrectDepthDependentOffset.expCorrRingSep;
                end
                axMeExp = obj.initaxis('Measured/expected');
                ringDistS1_good = ringDistS1(lGood);
                expCorrRingSep_good = expCorrRingSep(lGood);
                plot(axMeExp, ringDistS1_good(expCorrRingSep_good>=twoRingCutoff), expCorrRingSep_good(expCorrRingSep_good>=twoRingCutoff), ' ob')
                hold(axMeExp, 'on')
                plot(axMeExp, ringDistS1_good(expCorrRingSep_good<twoRingCutoff), expCorrRingSep_good(expCorrRingSep_good<twoRingCutoff), ' ok')
%                 plotSElink(axMeExp, ringDistS1_good(expCorrRingSep_good<43.7), expCorrRingSep_good(expCorrRingSep_good<43.7),ID_good(expCorrRingSep_good<43.7),se, ' ok')
                leftCor = max(axMeExp.XLim(1),axMeExp.YLim(1));
                rightCor =  max(axMeExp.XLim(2),axMeExp.YLim(2));
                plot(axMeExp, [leftCor rightCor], [leftCor rightCor],'-k')
                hold(axMeExp, 'off')
                title(axMeExp, 'Ring separation (nm)')
                xlabel(axMeExp, 'Measurement')
                ylabel(axMeExp, 'Expectation')
                
                if ~p.preview
                    warning('Z positions had been corrected before. The correction was not applied.')
                end
            else
                %% apply the correction
                if ~p.preview
                    obj.setPar('undoModule','CorrectDepthDependentOffset');
                    notify(obj.P,'backup4undo');
                    locs = obj.locData.loc;
                    if ~isfield(obj.locData.loc,'znm_rs_original')
                        obj.locData.loc.znm_rs_original = locs.znm;
                    end
    %                 for k = 1:length(sites_good)
    %                    sites_good(k).evaluation.SMLMModelFitGUI.allParsArg.value(idxRingDist) = newRingDist_good(k);
    %                    sites_good(k).evaluation.SMLMModelFitGUI.allParsArg.value(idxXrot) = rotAng(k,1);
    %                    sites_good(k).evaluation.SMLMModelFitGUI.allParsArg.value(idxYrot) = rotAng(k,2);
    %                    sites_good(k).evaluation.SMLMModelFitGUI.allParsArg.value(idxZ) = newZLower_good(k);
    %                 end

                    for k = 1:length(usedSites)
                        usedSites(k).evaluation.CorrectDepthDependentOffset.expCorrRingSep = newRingDist(k);
                        usedSites(k).evaluation.CorrectDepthDependentOffset.znm = z(k);
                    end
                    % Reset the filtering by znm
                    g = obj.getPar('mainGui');
                    znmFilter = [str2num(g.children.guiRender.children.Layer1.guihandles.znm_min.String) str2num(g.children.guiRender.children.Layer1.guihandles.znm_max.String)];
%                     d_znmFilter = relMeasuredOffsetFun(znmFilter);
                    znmFilter = correctedZFun(znmFilter);
                    g.children.guiRender.children.Layer1.guihandles.znm_min.String = num2str(znmFilter(1));
                    g.children.guiRender.children.Layer1.guihandles.znm_max.String = num2str(znmFilter(2));

                    obj.locData.loc.znm = correctedZFun(locs.znm);
%                     obj.locData.loc.znm = locs.znm-relMeasuredOffsetFun(locs.znm);           
                    obj.locData.regroup;
                    obj.locData.filter;
                    
                    for k = 1:length(usedSites)
                        usedSites(k).annotation.list3.value = 1;
                    end
                    
                    for k = find(~lGood)'
                        usedSites(k).annotation.list3.value = 7;
                    end
                end
            end
            
            
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function val = rotConversion(x,v)
    a = x(1);
    b = x(2);
    [x,y,z] = rotcoord3(0,0,-1, deg2rad(a), deg2rad(b), 0, 'XYZ');
    val = sqrt((v(1)-x)^2+(v(2)-y)^2+(v(3)-z)^2);
end

function pard=guidef

pard.preview.object=struct('String','Preview','Style','checkbox');
pard.preview.position=[2,1];
pard.preview.Width=1.5;

pard.plugininfo.type='ProcessorPlugin';


end