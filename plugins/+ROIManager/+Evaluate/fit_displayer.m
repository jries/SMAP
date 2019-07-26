classdef fit_displayer<interfaces.SEEvaluationProcessor

    methods
        function obj=fit_displayer(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        
        function out=run(obj,p)
            if p.reverse_z
                reverseFactor = -1;
            else
                reverseFactor = 1;
            end
            
            out=[];
            % make use of the SMLMModelFit class
            m1 = functionModel('C:\Users\ries\git\ries-private\SMLMModelFitter\models\CME3DContinuous.m');         % model 1
            fitting = SMLMModelFit({m1}, 'SolverName', 'fminsearchbnd','SolverOptions',{'Display','iter'});
            % fitting = SMLMModelFit({m1}, 'SolverName', 'particleswarm','SolverOptions',{'UseVectorized',true, 'Display','iter'});
            fitting.init                     % initiation
            fitting.roiSize = 500;
            fitResult = obj.site.evaluation.SMLMModelFit.fitResult;
            
            % update the arguments based on the fit results
            fitting.allParsArg = fitResult;
            
            % fetch locs
            locs = obj.getLocs({'xnm','ynm','znm','locprecnm','locprecznm'},'grouping','grouped','layer',1, 'size',[p.se_siteroi p.se_siteroi]);
            locs.xnm = locs.xnmrot;
            locs.ynm = locs.ynmrot;
            locs.znm = reverseFactor*locs.znm;
            
            % plot the result of DBSCAN clustering
            [~,re] = DBSCANsmap(locs, 30, 10);
            ind = re(:,4);
            indRange =1:1:max(ind)+1;
            [~,bestInd] = max(histcounts(ind(ind>0),indRange));
            indGood = ismember(ind,indRange(bestInd));
            indGood = indGood>0;
            ax1 = obj.setoutput('DBSCAN',1);
            plot3(ax1, locs.xnm(indGood), locs.ynm(indGood), locs.znm(indGood), ' or')
            hold(ax1,'on')
            plot3(ax1, locs.xnm(~indGood), locs.ynm(~indGood), locs.znm(~indGood), ' ob')
            xlabel(ax1, 'x')
            ylabel(ax1, 'y')
            zlabel(ax1, 'z')
            hold(ax1,'off')
            axis(ax1,'equal')
           
            % align the site based on the fit result
            lPars = fitting.exportPars(1,'lPar');
            locs2 = fitting.locsHandler(locs, lPars);
            
            %% Sections
            % get mPars
            mPars = fitting.exportPars(1,'mPar');
            closeAngle = mPars.closeAngle*pi/180;
            R = mPars.radius;           % radius of the sphere
            r = R*cos(closeAngle);      % radius of the open
            pixelSize = 2;
            theta = 0:pi/360:2*pi;      % full circle
            halfThickness=25;
            gaussFactor=0.4;
            
            % XY section
            [~,modP] = fitting.model{1}.modelFun(mPars,4); % offset determined by the model
            indXY = locs2.znm<halfThickness&locs2.znm>-halfThickness; % +-70 round the znm
            x = locs2.xnm(indXY);
            y = locs2.ynm(indXY);
            sigma = gaussFactor*locs2.locprecnm(indXY);
            img = renderImage(x,y,sigma,sigma);
            ax3 = obj.setoutput('Cross_section',1);
            imagesc(ax3, img)
            colormap(ax3, mymakelut('red hot'))
            
            
            
            % open
            xref = (r.*cos(theta)+fitting.roiSize/2)./pixelSize;
            yref = (r.*sin(theta)+fitting.roiSize/2)./pixelSize;
            
            % eqitorial
            xrefEq = (R.*cos(theta)+fitting.roiSize/2)./pixelSize;
            yrefEq = (R.*sin(theta)+fitting.roiSize/2)./pixelSize;
            
            hold(ax3,'on')
            plot(ax3, xref, yref, '- b')
            plot(ax3, xrefEq, yrefEq, '- g')
            hold(ax3,'off')
            
            subplot(2,2,1,ax3)
            axes(ax3)
            ax3_2 = subplot(2,2,2);
            axes(ax3)
            ax3_3 = subplot(2,2,3);
            
            % YZ section
            indYZ = locs2.xnm<halfThickness&locs2.xnm>-halfThickness; % +-70 round the xnm
            y = locs2.ynm(indYZ);
            z = locs2.znm(indYZ);
            sigma = gaussFactor*locs2.locprecnm(indYZ);
            sigmaZ = gaussFactor*locs2.locprecznm(indYZ);
            img = renderImage(y,z,sigma,sigmaZ);
            imagesc(ax3_2, img)
            colormap(ax3_2, mymakelut('red hot'))
            
            % close angle
            % open to the top
            thetaAnchor2 = pi-closeAngle;       % the reflection of the close angle in the 2nd quadrant
            theta = [closeAngle:-pi/360:-pi/2 3*pi/2:-pi/360:thetaAnchor2];
            theta = -theta;                     % reverse the circle (up-side-down)
            % after the reverse, now the open is to the bottom
            zrefClose = (R.*sin(theta)-modP.zOffset+fitting.roiSize/2)./pixelSize;
            yrefClose = (R.*cos(theta)+fitting.roiSize/2)./pixelSize;
            
            % open angle
            % open to the top
            theta = closeAngle:pi/180:thetaAnchor2;
            theta = -theta;                     % reverse the circle (up-side-down)
            % after the reverse, now the open is to the bottom
            zrefOpen = (R.*sin(theta) -modP.zOffset +fitting.roiSize/2)./pixelSize;
            yrefOpen = (R.*cos(theta)+fitting.roiSize/2)./pixelSize;
            
            hold(ax3_2,'on')
            plot(ax3_2, zrefClose, yrefClose, '- g')
            plot(ax3_2, zrefOpen, yrefOpen, '- b')
            hold(ax3_2,'off')
            set(ax3_2, 'Xdir','reverse')
            
            % ZX section
            indZX = locs2.ynm<halfThickness&locs2.ynm>-halfThickness; % +-70 round the ynm
            x = locs2.xnm(indZX);
            z = locs2.znm(indZX);
            sigma = gaussFactor*locs2.locprecnm(indZX);
            sigmaZ = gaussFactor*locs2.locprecznm(indZX);
            img = renderImage(z,x,sigmaZ,sigma);
            imagesc(ax3_3, img)
            colormap(ax3_3, mymakelut('red hot'))
            set(ax3_3, 'Ydir','normal')
            
            % the referece is just the axis-exchanged version as for the YZ
            % section
            hold(ax3_3,'on')
            plot(ax3_3, yrefClose, zrefClose, '- g')
            plot(ax3_3, yrefOpen, zrefOpen, '- b')
            hold(ax3_3,'off')
            
            axis(ax3,'equal')
            axis(ax3_2,'equal')
            axis(ax3_3,'equal')
            title(ax3, 'xy')
            title(ax3_2, 'yz')
            title(ax3_3, 'zx')
            
            %% point cloud
            ax2 = obj.setoutput('Point_cloud',1);
            xlabel(ax2, 'x')
            ylabel(ax2, 'y')
            zlabel(ax2, 'z')
            fitting.plot(locs, 'plotType', 'point', 'pixelSize', 2, 'axes',ax2);
            
            
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard = guidef(obj)
    pard.reverse_z.object = struct('Style','checkbox','String','Reverse Z', 'Value', 0);
    pard.reverse_z.position = [1 1]
    pard.reverse_z.Width = 1;
    
    pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
    pard.plugininfo.type='ROI_Evaluate';
end

function filterCallback(a,b,obj)
    if obj.guihandles.filtering.Value
        obj.guihandles.lockFilter.Enable = 'on';
    else
        obj.guihandles.lockFilter.Enable = 'off';
        obj.guihandles.lockFilter.Value = 0;
    end
end

function v = renderImage(a,b,sigmaA,sigmaB)
    roiSize = 500;
    bound = roiSize/2;
    rangex = [-bound bound]; % this should be linked to the ROI size
    rangey = rangex;
    rangez = rangex;
    pixelsize = 2;

    a=a(:)-rangey(1);b=b(:)-rangez(1);
    a=a/pixelsize; b=b/pixelsize;

    srec(1)=ceil((rangex(2)-rangex(1))/pixelsize);
    srec(2)=ceil((rangey(2)-rangey(1))/pixelsize);
    srec(3)=ceil((rangez(2)-rangez(1))/pixelsize);
    roiks=2.7
    G=creategausstemplate(roiks);

    [v,nlocs]=gaussrender_elliptc(single(a),single(b),uint32(srec),single(sigmaA),single(sigmaB),single(G.template),single(G.sigmatemplate),...
    single(roiks),single(ones(size(a))),int32(0),single(0),single(0), single([0 0]));
end

function gausstemplate=creategausstemplate(roiks) % create template
    % global gausstemplate
    % sigmatemplate=10;
    sizegauss=300;
    sigmatemplate=(sizegauss)/(2*roiks)/2; %for 2.5 sigma in both directions
    xg=-sizegauss:sizegauss;
    [Xg,Yg]=meshgrid(xg,xg);
    template=exp(-((Xg).^2+(Yg).^2)/2/sigmatemplate^2);
    gausstemplate.template=template;
    gausstemplate.sizegauss=sizegauss;
    gausstemplate.sigmatemplate=sigmatemplate;
end