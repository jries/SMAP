classdef continousChain
    properties
        name = {'r1', 'theta1', 'psi1',...
            'r2', 'theta2', 'psi2',...
            'r3', 'theta3', 'psi3',...
            'r4', 'theta4', 'psi4',...
            'r5', 'theta5', 'psi5'};
        fix = zeros([1 15]) ;
        value = repmat([30 0 0],[1 5]);
        lb = repmat([-15 -30 -30],[1 5]);
        ub = repmat([15 30 30],[1 5]);
        modelType = 'continuous';
        dimension = 3;
    end
    methods
        function obj = continousChain
        end
    end
    methods (Static)
        function [model, p]= reference(par, dx)
            
            model.x(1) = 0;
            model.y(1) = 0;
            model.z(1) = 0;
            
            for k = 1:numel(fieldnames(par))/3
                if k >1
                    par.(['theta' num2str(k)]) = par.(['theta' num2str(k)]) + par.(['theta' num2str(k-1)]);
                    par.(['psi' num2str(k)]) = par.(['psi' num2str(k)]) + par.(['psi' num2str(k-1)]);
                end
                delx = par.(['r' num2str(k)]).*cos(par.(['theta' num2str(k)])*180/pi);
                delz = par.(['r' num2str(k)]).*sin(par.(['theta' num2str(k)])*180/pi);
                [delx,dely] = rotcoord(delx,0,par.(['psi' num2str(k)])*180/pi);
                model.x(k+1) = model.x(k)+delx;
                model.y(k+1) = model.y(k)+dely;
                model.z(k+1) = model.z(k)+delz;
            end
            
            [pp,p]=csaps(model.x,[model.y;model.z]);
            model.x = linspace(min(model.x),max(model.x),range(model.x)/dx);
            val=fnval(pp,linspace(min(model.x),max(model.x),range(model.x)/dx));
            model.y = val(1,:);
            model.z = val(2,:);
            model.n = ones(size(model.x));
            p = [];
        end
    end
end