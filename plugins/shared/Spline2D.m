classdef Spline2D<handle
    properties
        max_i;
        coeff;
    end
    methods
        function obj = Spline2D(d)
%         function obj = Spline2D( d, coeff = false, verbose = false)
            d_size = size(d,1);
            obj.max_i = d_size-1;
            ys = [];
            for i = 1:d_size
                ys = [ys Spline1D(d(i,:))];%append(spline1D.Spline1D(d[i,:]))
            end
            xs = [];
            cx = 0;
            %For obtainning the 16 coefficients, we need 16 equations (16
            %known points); therefore, we resample each grid 16 times with
            %spacing 1/3 in x and y;
            while cx <= (obj.max_i+0.01)
                if cx > obj.max_i
                    cx = obj.max_i;
                end
                xv = zeros(d_size,1);
                for i = 1:d_size
                    xv(i) = ys(i).f(cx);
                end
                xs = [xs Spline1D(xv)];
                cx = cx + 1/3;
            end
            obj.coeff = zeros(obj.max_i, obj.max_i,16);
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
%             for i = 1:(d_size-1)
%                 for j = 1:(d_size-1)
            for i = 1:obj.max_i
                for j = 1:obj.max_i
                    for k = 1:4
                        sp = xs(3*(i-1) + k);
                        for l = 1:4
                            cy = j-1 + (l-1)/3;
                            b((k-1)*4+l)=sp.f(cy);
                        end
                    end
%                     b1=b(1)
%                     d1=d(j,i)
%                     b2=b(4)
%                     d2=d(j+1,i)
%                     b3=b(13)
%                     d3=d(j,i+1)
%                     b4=b(16)
                    
                    
                   
%                     d4=d(j+1,i+1)
                    x = A\b;
                    obj.coeff(i,j,:) = reshape(x(:),[1,1,16]);

                end
            end            
        end
        

        function yval = f(obj,y,x)
%             y=1;
%             x=1;
  % The first pixel is 0 regarding the function roundAndCheck
            [ix, x_diff] = roundAndCheck(x,obj.max_i);
            [iy, y_diff] = roundAndCheck(y,obj.max_i);
            if ix == -1|| iy == -1
                yval = 0;
                return
            end
            yval = 0;
            for i = 1:4
                for j = 1:4
                    yval = yval+obj.coeff(ix,iy,4*(i-1)+j)*x_diff^(i-1)*y_diff^(j-1);
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