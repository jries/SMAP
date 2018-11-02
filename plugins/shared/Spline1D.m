classdef Spline1D<handle
    properties
        max_i;
        coeff;
    end
    methods
        function obj = Spline1D(y)
            if nargin ==0
                return
            end
            obj.max_i = numel(y)-1;
            
            b = zeros(numel(y),1);
            %??
            for i = 1:(numel(y)-2)
                b(i+1) = 6*(y(i)-2*y(i+1)+y(i+2));
            end
            A = zeros(numel(y),numel(y));
            A(1,1)=1;
            A(end,end)=1;
            for i = 1:(numel(y)-2)
                if 0
                else
                    A(i+1,i)=1;
                    A(i+1,i+1)=4;
                    A(i+1,i+2)=1;
                end
            end
            %??
            M = A\b;
            
            %compute spline coefficients
            % f = @(x) a(i) + bj(i).*(x-xi(i)) + cj(i).*(x-xi(i)).^2 + dj(i).*(x-xi(i)).^3;
            obj.coeff = zeros(numel(y)-1,4);
            for i = 1:(numel(y)-1)
                obj.coeff(i,4) = (M(i+1)-M(i))/6;%d(i)
                obj.coeff(i,3) = M(i)/2;%c(i)
                obj.coeff(i,2) = (y(i+1)-y(i))-(M(i+1)+2*M(i))/6;%b(i)
                obj.coeff(i,1) = y(i);%a(i)
                
            end
            
        end
        
        function [ix, x_diff] = roundAndCheck(obj, x, max_x)
            if x<0||x>max_x
                disp(['value out of range'])
                ix = -1;
                x_diff = -1;
                return
            end
            x_floor = floor(x);
            x_diff = x - x_floor;
            ix = floor(x_floor)+1;
            if x == max_x
                ix = ix-1;
                x_diff = 1;
            end
        end
        
        function yval = f(obj,x)
             max_x=obj.max_i;
             if x<0||x>max_x
                disp(['value out of range'])
                ix = -1;
                x_diff = -1;
%                 return
             else
            x_floor = floor(x);
            x_diff = x - x_floor;
            ix = floor(x_floor)+1;
            if x == max_x
                ix = ix-1;
                x_diff = 1;
            end
             end
            
%             [ix, x_diff] = obj.roundAndCheck(x,obj.max_i);
            if ix == -1
                yval = 0;
                return
            end
            yval = 0;
            coeff=obj.coeff;
            for i = 1:4
                yval = yval + coeff(ix,i)*x_diff^(i-1);
            end
            
        end
        
        
        
    end
    
end























