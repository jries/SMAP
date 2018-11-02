classdef Spline2D_v3<handle
    properties
        max_ix;
        max_iy;
        coeff;
    end
    methods
        function obj = Spline2D_v3(d)
            if nargin ==0
                return
            end
%         function obj = Spline2D( d, coeff = false, verbose = false)
            d_sizex = size(d,1);
            d_sizey = size(d,2);
            obj.max_ix = d_sizex-1;
            obj.max_iy = d_sizey-1;
%             ys = [];
            for i = d_sizex:-1:1
%                 ys = [ys Spline1D(d(i,:))];%append(spline1D.Spline1D(d[i,:]))
                 ys(i) = Spline1D(d(i,:));%append(spline1D.Spline1D(d[i,:]))
            end
            xs = [];
            cx = 0;
            %For obtainning the 16 coefficients, we need 16 equations (16
            %known points); therefore, we resample each grid 16 times with
            %spacing 1/3 in x and y;
            while cx <= (obj.max_iy+0.01)
                if cx > obj.max_iy
                    cx = obj.max_iy;
                end
                xv = zeros(d_sizex,1);
                for i = 1:d_sizex
                    xv(i) = ys(i).f(cx);
                end
                xs = [xs Spline1D(xv)];
                cx = cx + 1/3;
            end
            obj.coeff = zeros(obj.max_ix, obj.max_iy,16);
            A = zeros(16,16);
            for i = 1:4
                dx = (i-1)/3;
                for j = 1:4
                    dy = (j-1)/3;
                    for k = 1:4
                        for l = 1:4
                            A((i-1)*4+j,(k-1)*4+l) = dx^(k-1) * dy^(l-1);
                        end
                    end
                end
            end
            
            b = zeros(16,1);

%             for i = 1:obj.max_iy
%                 for j = 1:obj.max_ix
%                     for k = 1:4
%                         sp = xs(3*(i-1) + k);
%                         for l = 1:4
%                             cy = j-1 + (l-1)/3;
%                             b((k-1)*4+l)=sp.f(cy);
%                         end
%                     end
% 
%                     x = A\b;
%                     %need check
%                     obj.coeff(j,i,:) = reshape(x(:),[1,1,16]);
% 
%                 end
%             end
            
            for i = 1:obj.max_ix
                for j = 1:obj.max_iy
                    for k = 1:4
                        sp = xs(3*(j-1) + k);
                        for l = 1:4
                            cy = i-1 + (l-1)/3;
                            b((k-1)*4+l)=sp.f(cy);
                        end
                    end
                    x = A\b;
                    %need check
                    obj.coeff(i,j,:) = reshape(x(:),[1,1,16]);
                    
                end
            end
            
            
        end
        

        function yval = f(obj,x,y)
%             y=1;
%             x=1;
  % The first pixel is 0 regarding the function roundAndCheck
            [ix, x_diff] = roundAndCheck(x,obj.max_ix);
            [iy, y_diff] = roundAndCheck(y,obj.max_iy);
            if ix == -1|| iy == -1
                yval = 0;
                return
            end
            yval = 0;
            coeff=obj.coeff;
            for i = 1:4
                for j = 1:4
                    yval = yval+coeff(ix,iy,4*(i-1)+j)*y_diff^(i-1)*x_diff^(j-1);
                end
            end
        end
        
        
        
    end
    
end

function [ix, x_diff] = roundAndCheck(x, max_x)
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