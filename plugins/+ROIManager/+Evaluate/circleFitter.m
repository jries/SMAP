classdef circleFitter<interfaces.SEEvaluationProcessor
    properties
        boundary
    end
    methods
        function obj=circleFitter(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, inp)
            % Import the kimograph
%             site = obj(1).locData
            layers=find(inp.sr_layerson);
            locs=obj.getLocs({'locprecnm','PSFxnm','xnm','ynm','frame'},'layer',1,'size',inp.se_siteroi);  
            figure(30)
            scatter(locs.xnm, locs.ynm);
            test(test>=1)=1;
            
            
            [x, y] = meshgrid(400:800, 400:800);
            scatter3(x(:), y(:), repelem(0, 40401)')
            
            [xq, yq] = meshgrid(locs.xnm, locs.ynm);
            
            zq = griddata(x(:), y(:), repelem(0, 401*401)', xq ,yq);
            
            contour(xq ,yq, zq)
            
            slideStep = inp.slideStep;
            gridRange = inp.gridRange;
            adjM = inp.adjM;

            % expand test
            oriSize = size(test);
            expand = zeros(oriSize(1),gridRange*slideStep);
            test = [test, expand];

            % Binarize the kimograph
            [testRefx, testRefy] = find(test);
            [testNZx, testNZy] = find(test);

            % use euclidean distances to filter out noises
            K = round(length(testRefx)*5/100);

            ED = pdist2([testNZx testNZy],[testNZx testNZy],'euclidean');
            sortedED = sort(ED,2);
            sortedED = sortedED(:,1:K);
            rmean = mean(sortedED,2);

            EDcutoff = quantile(rmean, 0.95);

            test(sub2ind(size(test), testNZx([rmean>=EDcutoff]), testNZy([rmean>=EDcutoff]))) = 0;

            % Get cordinates of the kemograph
            Size = size(test);
            [CorX CorY] = meshgrid(1:Size(1), 1:Size(2));
            Order = sub2ind(Size, CorX(:), CorY(:));

            % Make grids (see gridRange)
            adjCorX = ceil(CorX(:)./gridRange);
            adjCorY = ceil(CorY(:)./gridRange);

            % Get density
            D = accumarray([adjCorX adjCorY], test(Order), [], @sum);

            % Filter out grids that might be background
            cutoff = 3; %#Par
            D(D<=cutoff)=0;

            % Do sliding subtraction (see slideStep)
            % to see the next slideStep-1 right neighbor is higher or not than the
            % target
            SizeD = size(D);

            diff = cell(slideStep-1,1);
            diffAll = zeros(SizeD(1), SizeD(2)-slideStep);
            for i = 1:(slideStep-1)
                diffAll = diffAll + (D(:,1+i:end-slideStep+i) - D(:,1:(end-slideStep)) < 0);
            end

            % Get grids with 4 zero grids on its right-hand side
            diffAllLog = diffAll == slideStep-1; 
            [bX,bY] = find(diffAllLog, 1, 'last');

            % Get boundary (itself is 1, and its cumsum is also 1)
            boundary = diffAllLog & cumsum(diffAllLog,2,'reverse') == 1;

            %figure(50)
            %mesh(boundary)
            %rotate3d on

            % Find the corrdinates of boundary
            [Cx, Cy]= find(boundary,sum(boundary(:)),'last');

            % Make sure every frame (x-axies) have a point
            bounSize = size(boundary);
            CxFinal = (1:bounSize(1))';
            CyFinal = zeros(bounSize(1),1);
            CyFinal(Cx) = Cy;

            % Find the cummax (next points should always >= the previous points)
            CyFinal = cummax(CyFinal);
            CxGridParent = CxFinal*gridRange;
            CyGridParent = CyFinal*gridRange;

            meaningfulX = (CxGridParent+1)<=Size(1);
            CxGridParent = CxGridParent(meaningfulX);
            CyGridParent = CyGridParent(meaningfulX);

            CxFParent = (0:Size(1)-1)';
            CyFParent = zeros(Size(1),1);
            CyFParent(CxGridParent+1) = CyGridParent;
            CyFParent = cummax(CyFParent);


            %figure(50)
            %scatter(testNZx, testNZy)
            %hold on
            %plot(CxFParent, CyFParent)


            Mref = calMeasurement(CxFParent,CyFParent,testNZx,testNZy,Size);
            MrefO = Mref;
            CyFParentA = CyFParent;
            for i=size(CxFParent):-1:1
               thisFrameY = testNZy(testNZx == CxFParent(i));
               thisFrameY = thisFrameY(thisFrameY>CyFParent(i));
               ii = 1;
               futherStep = 0;
               while ii <= length(thisFrameY)
                   CyFParentO = CyFParentA;

                   % New boundary
                   CyFParentA(i)=thisFrameY(ii);

                   CyFParentA = cummax(CyFParentA);

                   % fit to new boundary
                   Mnew = calMeasurement(CxFParent,CyFParentA,testNZx,testNZy,Size);
                   if Mnew >= MrefO*adjM && Mnew >= Mref %#Par
                       MrefO = Mnew;
                       futherStep = 0;
                   else
                       if futherStep <= 5
                        CyFParentA = CyFParentO;
                        Mnew = MrefO;
                        futherStep = futherStep+1;
                       else
                        CyFParentA = CyFParentO;
                        Mnew = MrefO;
                        break
                       end
                   end
                   ii=ii+1;
               end
                if isinteger(i/10)
                    figure(51)
                    scatter(testNZx, testNZy)
                    hold on
                    plot(CxFParent, CyFParentA)
                    hold off
                end
            end
            
            KimoPar = obj.site.evaluation.BALM_fibril_growth;
            fig = KimoPar.kimograph;
            co=quantile(fig(:),0.999);
            fig(fig>co)=co;
            h=obj.setoutput('kimograph');
            imagesc(h,(fig))
            %xlabel(h,'xnm')
            %ylabel(h,'frame')
            %xticklabels(h, KimoPar.xn)
            %yticklabels(h, KimoPar.fr)
            hold(h,'on')
            plot(h,CyFParentA,CxFParent, 'LineWidth', 1, 'Color', 'w')
            hold(h,'off')

            out.boundary = [CxFParent CyFParentA];
        end
     
        function pard=guidef(obj)
            pard=guidef;
        end
    end

end



function pard=guidef
pard.t_gridRange.object=struct('Style','text','String','Bin size of the grids');
pard.t_gridRange.position=[1,1];
pard.t_gridRange.Width=2;

pard.gridRange.object=struct('Style','edit','String',5);
pard.gridRange.position=[1,3];
pard.gridRange.TooltipString = 'If you set it as 5, it means before the density comparison, every grid will be set to cover a 5-by-5 area in the original coordinates';
            
pard.t_slideStep.object=struct('Style','text','String','Slide step(s)');
pard.t_slideStep.position=[2,1];
pard.t_slideStep.Width=2;

pard.slideStep.object=struct('Style','edit','String',5);
pard.slideStep.position=[2,3];
pard.slideStep.TooltipString = 'If you set it as 5, it means during the density comparison, every grid value will be campared to its following 4 (5 minus 1, which means the reference grid itself) right neighbors';

pard.t_adjM.object=struct('Style','text','String','Adjustment of M');
pard.t_adjM.position=[3,1];
pard.t_adjM.Width=2;

pard.adjM.object=struct('Style','edit','String',1.0003);
pard.adjM.position=[3,3];
pard.adjM.TooltipString = 'If you set it as 1.0003, it means during the optimization, if the measurment of current step (Mcur) is 0.0003-time worse than the measurment of the previous step (Mpre), Mcur will still be considered as a good result. The measurment, which defines the boundary is good or not, can be definde by users.';

% pard.dxt.Width=3;
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
pard.plugininfo.type='ROI_Evaluate';
end

function M = calMeasurement(x,y,qx,qy,Size)
    ref = interp1(x,y,qx);
    rightIdx = qy > ref;
    leftIdx = ~rightIdx;
    leftA = sum(y);
    rightA = Size(1)*Size(2)-sum(y);
    leftD = sum(leftIdx)/leftA;
    rightD = sum(rightIdx)/rightA;
    M = 0-((leftD-1)^2 + (rightD-0)^2)^(1/2);
end