classdef CME2Dfitting2<interfaces.SEEvaluationProcessor
    methods
        function obj=CME2Dfitting2(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function makeGui(obj)
            makeGui@interfaces.SEEvaluationProcessor(obj);
%             obj.guihandles.saveimagesb.Callback={@saveimagesb_callback,obj};
        end
        function out=run(obj,p)
            %% Parameters
            if p.method.Value ==2               %p.method.Value is the corresponding value of the chosen method.
                global capmatrix bottommatrix
            end
            
            if p.method.Value==1
                p.onStep = [];
            else
                p.onStep = [p.method.Value];
            end
            
            if p.plot
                p.onStep = [p.onStep 99];
            end
            
            out = [];
            
            %% Initiation
            if 0
                imgSize = 301;
                [X,Y] = meshgrid(-(imgSize-1)/2:(imgSize-1)/2,-(imgSize-1)/2:(imgSize-1)/2);
                test = timimgStructure([ 0, 100, -50, 1, 20,70, 20,40, 30],X(:),Y(:));

                [X,Y] = meshgrid(1:imgSize,1:imgSize);
                img = zeros(imgSize);
                ind = sub2ind([imgSize imgSize], Y(:), X(:));
                img(ind) = test;
                figure; imagesc(img)
            end
            
            info = [];
            locs = obj.getLocs({'xnmrot','ynmrot', 'channel'},'size',obj.P.par.se_siteroi.content/2,'grouping','grouped', 'layer', 1);
            if p.filterByMeanshift
                locs = structfun(@(x)x(obj.site.evaluation.CME2CSide_pg.filter'),locs, 'UniformOutput', 0);
            end
            pxSize = 5;
            binRange = -250:pxSize:250;
            
            %% Main part
            for k = 1:length(p.onStep)
                switch p.onStep(k)
                    case 2                          % Original implementation
                        [fittedVals, logLikli]= MEMLETCL([locs.xnmrot locs.ynmrot], 'nucleationNCoat(x, y, xcenter, ycenter, distance, ampAll, ampCap, angle, innerRadR, outerRadR, innerRadC, outerRadC, thickness, pxsize, roisize)', 'x,y', 'xcenter,ycenter,distance,ampAll,ampCap,angle,innerRadR,outerRadR,innerRadC,outerRadC,thickness,pxsize,roisize', '-20,30,-250,1,25,-45,30,70,20,45,60,5,300','20,100,0,1,400,45,30,70,20,45,60,5,300', ['0,-20,', num2str(-obj.site.evaluation.CME2CSide_pg.estimatedZRangeByDensity/2), ',1,100,0,30,70,20,45,60,5,500'],25);
                        c = fittedVals(1:6)';
                    case 3                          % New single-colour mode: move the coordinates instead of the model
                        capmatrix = obj.getPar('capmatrix');
                        bottommatrix = obj.getPar('bottommatrix');
                        [fittedVals, logLikli]= MEMLETCL([locs.xnmrot locs.ynmrot zeros(size(locs.ynmrot))], 'nucleationNCoat2(x, y, z, xShift, distanceTop, distanceBottom, ampAll, ampCap, angle, innerRadR, outerRadR, innerRadC, outerRadC, thickness, pxsize, roisize, varargin)', 'x,y,z', 'xShift,distanceTop,distanceBottom,ampAll,ampCap,angle,innerRadR,outerRadR,innerRadC,outerRadC,thickness,pxsize,roisize', '-20,20,-150,1,100,-45,30,70,20,45,60,5,300','20,150,0,1,100,45,30,70,20,45,60,5,300', ['0,-20,', num2str(-obj.site.evaluation.CME2CSide_pg.estimatedZRangeByDensity/2), ',1,100,0,30,70,20,45,60,5,500'],25);
                        c = fittedVals(1:6)';
                    case 4                          % Particle swarm optimization
                        if ~p.dualColour                % Switch between single-/dual-colour mode 
                            ch = 0;
                            info.mode = '1C';
                        else
                            ch = locs.channel;
                            info.mode = '2C';
                        end
                        
                        objfun =@(x) -sum(log10(nucleationNCoat3(locs.xnmrot, locs.ynmrot, ch, x(1), x(2), x(3), x(4), x(5),x(6),x(7),x(8),x(9),x(10), 'SumOffset', p.sumOffset, 'RoiSize', 300)));
%                         options = optimoptions('surrogateopt','PlotFcn','surrogateoptplot',...
%                                  'InitialPoints',[0 20 num2str(-obj.site.evaluation.CME2CSide_pg.estimatedZRangeByDensity/2) 0],'MaxFunctionEvaluations',5000);
                        fittedVals = particleswarm(objfun,10, p.lb, p.ub);
                        
                        [~,a1,a2,ampC, nLocs]= nucleationNCoat3(locs.xnmrot, locs.ynmrot, ch, fittedVals(1), fittedVals(2), fittedVals(3), fittedVals(4),fittedVals(5),fittedVals(6),fittedVals(7),fittedVals(8),fittedVals(9),fittedVals(10), 'SumOffset', p.sumOffset, 'RoiSize', 300);
                        c = [fittedVals(1:3) 1 fittedVals(4:10)];
                        info.method = 'PSO_likelihood';
                    case 5                          % Particle swarm optimization
                        if ~p.dualColour                % Switch between single-/dual-colour mode 
                            ch = 0;
                            info.mode = '1C';
                        else
                            ch = locs.channel;
                            info.mode = '2C';
                        end
                        
                        objfun =@(x) -sum(nucleationNCoat3(locs.xnmrot, locs.ynmrot, ch, x(1), x(2), x(3), x(4), x(5),x(6),x(7),x(8),x(9),x(10),'SumOffset', p.sumOffset, 'roisize', 300));
%                         options = optimoptions('surrogateopt','PlotFcn','surrogateoptplot',...
%                                  'InitialPoints',[0 20 num2str(-obj.site.evaluation.CME2CSide_pg.estimatedZRangeByDensity/2) 0],'MaxFunctionEvaluations',5000);
                        fittedVals = particleswarm(objfun,10, p.lb, p.ub);
                        
                        [~,a1,a2]= nucleationNCoat3(locs.xnmrot, locs.ynmrot, ch, fittedVals(1), fittedVals(2), fittedVals(3), fittedVals(4),fittedVals(5),fittedVals(6),fittedVals(7),fittedVals(8),fittedVals(9),fittedVals(10));
                        c = [fittedVals(1:3) 1 fittedVals(4:10)];
                        info.method = 'PSO_xcorr';
                        
                        
                        
                    case {6, 7}                          % Particle swarm optimization
                        if ~p.dualColour                % Switch between single-/dual-colour mode 
                            ch = 0;
                            info.mode = '1C';
                        else
                            ch = locs.channel;
                            info.mode = '2C';
                        end
                        
                        capmatrix = obj.getPar('capmatrix');
                        bottommatrix = obj.getPar('bottommatrix');
                        xcor =  obj.getPar('xcor');
                        ycor =  obj.getPar('ycor');
                        if isempty(capmatrix)
                            [capmatrix, bottommatrix] = generatePDFImg([0 0 30 70 20 40 60], 300, 10, 5);
                            [xcor,ycor]=meshgrid(1:length(capmatrix), 1:length(capmatrix));
                            obj.setPar("capmatrix", capmatrix);
                            obj.setPar("bottommatrix", bottommatrix);
                            obj.setPar("xcor", xcor);
                            obj.setPar("ycor", ycor);
                        end
                        
                        objfun =@(x) -sum(log10(nucleationNCoat2_2(locs.xnmrot, locs.ynmrot, ch, x(1), x(2), x(3), x(4), x(5),x(6),x(7),x(8),x(9),x(10),'gridX', xcor,'gridY',ycor, 'capmatrix',capmatrix,'bottommatrix',bottommatrix, 'roisize', 300, 'SumOffset', p.sumOffset)));
%                         options = optimoptions('surrogateopt','PlotFcn','surrogateoptplot',...
%                                  'InitialPoints',[0 20 num2str(-obj.site.evaluation.CME2CSide_pg.estimatedZRangeByDensity/2) 0],'MaxFunctionEvaluations',5000);
                        
                    switch p.onStep(k)
                        case 6
                            fittedVals = particleswarm(objfun,10, p.lb, p.ub);
                        case 7
                            options = optimoptions('surrogateopt','Display','off', 'MaxFunctionEvaluations',300, 'PlotFcn', []);
                            fittedVals = surrogateopt(objfun, p.lb, p.ub, options);
                    end
                        [~,a1,a2,ampC, nLocs]= nucleationNCoat2_2(locs.xnmrot, locs.ynmrot, ch, fittedVals(1), fittedVals(2), fittedVals(3), fittedVals(4),fittedVals(5),fittedVals(6),fittedVals(7),fittedVals(8),fittedVals(9),fittedVals(10),'gridX', xcor,'gridY',ycor, 'capmatrix',capmatrix,'bottommatrix',bottommatrix, 'roisize', 300, 'SumOffset', p.sumOffset);
                        c = [fittedVals(1:3) 1 fittedVals(4:10)];
                        info.method = 'gaussianFiltered';
                        

                    case 99                         % Visulization
                        if ~exist('c')
                            c = obj.locData.SE.currentsite.evaluation.CME2Dfitting2.parameters(1:11);
                            a1 = obj.locData.SE.currentsite.evaluation.CME2Dfitting2.parameters(12);
                            a2 = obj.locData.SE.currentsite.evaluation.CME2Dfitting2.parameters(13);
                            ampC = obj.locData.SE.currentsite.evaluation.CME2Dfitting2.parameters(14);
                            nLocs = obj.locData.SE.currentsite.evaluation.CME2Dfitting2.parameters(15);
                        end
                        if ismember(2, p.onStep)
                            c = transPar(c);
                        end
                        

        %                 figure(110); imagesc(imgfit);
        %                 figure(111); imagesc(img);
        %                 figure(112); imagesc(img-imgfit./sum(imgfit(:))*length(locs.xnmrot));
                        if c(5) ==0
                            c(5) = ampC;
                        end
                        inp.px = 1;
                        [gridx, gridy, gridv] = gridTimimgStructure(c,'PixelSize',inp.px,'RoiSize',300);
                        img = histcounts2(locs.ynmrot, locs.xnmrot, -150:inp.px:150, -150:inp.px:150);
                        h = fspecial('gaussian',round(30/inp.px),round((30/inp.px)/2));
                        h2 = fspecial('gaussian',round(40/inp.px),round((40/inp.px)/2));
                        img = filter2(h,img);
                        imgfit = filter2(h2,gridv);
                        fig1 = figure(113);
                        ax1 = axes(fig1);
                        cla(ax1);        
                        imagesc(ax1, imgfit);
                        pbaspect(ax1, [1 1 1]);
                        colormap(ax1,mymakelut('red'))
                        l1 = getframe(ax1);
                        l1 = l1.cdata;
                        
                        fig4 = figure(116);
                        ax4 = axes(fig4);
                        cla(ax4);        
%                         imagesc(ax4, imgfit,[0 2.2e-10]);
                        imagesc(ax4, imgfit);
                        pbaspect(ax4, [1 1 1]);
                        colormap(ax4,mymakelut('gray'))

                        fig2 = figure(114);
                        ax2 = axes(fig2);
                        cla(ax2);
                        imagesc(ax2, img);
                        pbaspect(ax2, [1 1 1]);
                        colormap(ax2,mymakelut('green'))
                        l2 = getframe(ax2);
                        l2 = l2.cdata;

                        imgall = l1+l2;
                        fig3 = figure(115);
                        ax3 = axes(fig3);
                        cla(ax3);                        
                        imagesc(imgall)
                        pbaspect(ax3, [1 1 1]);
                end
                
                out.parameters=[c a1 a2 ampC nLocs];
                out.info = info;
            end
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
pard.t1.object = struct('String', 'Method', 'Style','text');
pard.t1.position =[1,1]; 
pard.t1.Width = 1;

pard.method.object=struct('String',{{'None','Single_old','Single_rotCoord','PSO_likelihood','PSO_xcorr', 'gaussianFiltered_PSO', 'gaussianFiltered_Sur'}},'Style','popupmenu', 'Value', 1);
pard.method.position=[1,2];
pard.method.Width=1;

pard.plot.object=struct('String','Plot','Style','checkbox','Value',1);
pard.plot.position=[2,1];
pard.plot.Width=1;

pard.dualColour.object=struct('String','Dual-colour','Style','checkbox','Value',1);
pard.dualColour.position=[2,2];
pard.dualColour.Width=1;

pard.t2.object = struct('String', 'Lower bound', 'Style','text');
pard.t2.position =[3,1]; 
pard.t2.Width = 1;

pard.lb.object = struct('String', '-60 30 -200 0.33 -45 30 70 20 40 60', 'Style','edit');
pard.lb.position =[3,2]; 
pard.lb.Width = 2;

pard.t3.object = struct('String', 'Upper bound', 'Style','text');
pard.t3.position =[4,1]; 
pard.t3.Width = 1;

pard.ub.object = struct('String', '60 120 0 3 45 30 70 20 40 60', 'Style','edit');
pard.ub.position =[4,2]; 
pard.ub.Width = 2;

pard.t2.object = struct('String', 'Offset', 'Style','text');
pard.t2.position =[5,1]; 
pard.t2.Width = 2;

pard.sumOffset.object = struct('String', '0.1', 'Style','edit');
pard.sumOffset.position =[5,2]; 
pard.sumOffset.Width = 2;

pard.filterByMeanshift.object = struct('String', 'Filter(mean-shift)', 'Style','checkbox', 'Value', 1);
pard.filterByMeanshift.position =[6,1]; 
pard.filterByMeanshift.Width = 1.5;

pard.plugininfo.type='ROI_Evaluate';
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};

end
function newC = transPar(c)
    newC = [-c(1) -c(3) -(c(2)-c(3)) c(4) c(5) c(6)];
end

function [imgC, imgR]= generatePDFImg(c,roisize, size, std)
% default c: [0 0 x(6) x(7) x(8) x(9) x(10)]
imgC = zeros([roisize+1 roisize+1]);
imgR = zeros([roisize+1 roisize+1]);

[x, y] = meshgrid(1:roisize+1, 1:1:roisize+1);
indImg = sub2ind([roisize+1 roisize+1], x(:), y(:));

v1 = thickRing([c(1)+roisize/2+1, c(2)+roisize/2+1, c(3), c(4), c(7)], x(:), y(:));
% v1 = 4([c(1)+roisize/2+1, c(2)+roisize/2+1, c(3), c(4), c(7)], x(:), y(:));
v2 = cap([c(1)+roisize/2+1, c(2)+roisize/2+1, c(5), c(6), c(7)], x(:), y(:));

imgR(indImg) = v1;
imgC(indImg) = v2;
h = fspecial('gaussian', size, std);
imgR = filter2(h,imgR)/sum(imgR(:));
imgC = filter2(h,imgC)/sum(imgC(:));
end