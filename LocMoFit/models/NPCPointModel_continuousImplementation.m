classdef NPCPointModel_continuousImplementation
    properties
        name = {'ringDistance', 'azimuthalShift', 'diameter'};
        fix = [0 0 1] ;
        value = [0 0 53.7];
        lb = [0 -inf 0];
        ub = [70 inf 0];
        min = [0 -inf 0];
        max = [70 inf 100];
        modelType = 'continuous'
        dimension = 3;
    end
    methods
        function obj = NPCPointModel
        end
    end
    methods (Static)
        function [model, p]= reference(par, dx)
        % model: transformed locs
        % p: parameters for the model

        numParticle = max(structfun(@length, par));
        % need to open the parameters as options

        % radius=50;
        % lpc=4; %number of corners
        % dth=pi/48; %shift between corners
        % theta=0:pi/4:2*pi-pi/1000; %8 corners
        % thetaa=[];
        % for k=1:lpc
        %     thetaa=[thetaa theta+k*dth];
        % end
        % [model.x,model.y]=pol2cart(thetaa,radius);
        % model.channel=ones(size(model.x));
            zwobble=0;
            rwobble=0;

        % ring_angle=13; %deg %Nup107
            ring_angle=par.azimuthalShift; 
        % ring_angle=8; %deg %Nup96
            pair_angle=12;
            dz=par.ringDistance+randn.*zwobble;
            corners=8;
            r = zeros(numParticle, 4);
            r(:,1)=par.diameter';                          % first copy
            r(:,2)=par.diameter';                          % second copy
        % r(1)=53.7;                              % first copy
        % r(2)=54;                                % second copy
            r(:,4)=r(:,1);r(:,3)=r(:,2);
            r=r+randn.*rwobble;
            a = zeros(numParticle, 4);
            a(:,1)=0;                                 % the first copy (ring A)
            a(:,3)=a(:,1)+pair_angle';                   % the second copy (ring A)
            a(:,2)=a(:,1)+ring_angle';                   % the first copy (ring B)
            a(:,4)=a(:,2)+pair_angle';                   % the second copy (ring B)
            z = zeros(numParticle, 4);
            z(:,1)=dz'/2;
            z(:,2)=-dz'/2;z(:,4)=-dz'/2;z(:,3)=dz'/2;
        %NUP 96
        % rN(1)=43.2;
        % rC(1)=47.2;
        % rN(2)=49.1;
        % rC(2)=52.5;

        %total angle between corresponding points as in excel sheet, deg

        % aN(1)=9.1;
        % aC(1)=15.9;
        % aN(2)=15;
        % aC(2)=7.8;

        % r=rC;
        % a=aC;
        %to deg, angle from ref:
        % a=a/2/180*pi;
            a=a/180*pi;


        % a(3)=-a(1);a(4)=-a(2);
        % z(1)=0;z(2)=50;z(3)=50;z(4)=0;

            rall=zeros(size(r,2),corners, size(r,1));
            angall=zeros(size(r,2),corners, size(r,1));
            zall=zeros(size(r,2),corners, size(r,1));

            dphi=0:2*pi/corners:2*pi-pi/1000;

            for k=1:size(r,2)
                rall(k,:,:)=(r(:,k)+repelem(0,corners))';
                angall(k,:,:)=(a(:,k)+dphi)';
                zall(k,:,:)=(z(:,k)+repelem(0,corners))';
            end

            model.x = zeros(size(rall,1)*size(rall,2),size(rall,3));
            model.y = model.x;
            model.z = model.x;
            model.channel = model.x;
            for k = 1:size(rall,3)
                tempAng = angall(:,:,k);
                tempR = rall(:,:,k);
                [model.x(:,k),model.y(:,k)]=pol2cart(tempAng(:),tempR(:));
                tempZ = zall(:,:,k);
                model.z(:,k)=tempZ(:);
                model.channel(:,k)=ones(size(model.x,1),1);
            end
            model.n = model.channel; % should be changed properly later

        %     likelihood = gaussDist(model,double(locs.ynmrot), double(locs.xnmrot), double(locs.znm), 8,8,10);
        %     likelihoodBoun = gaussDist(model,0,0,0, 16,16,20);

            p.ring_angle=ring_angle;
            p.pair_angle=pair_angle;
            p.ring_distance=dz;
            p.radius=r;
        end
    end
end
