classdef Spline3D<handle
    properties
        max_i;
        coeff;
    end
    methods
        function obj = Spline3D(d)
            if size(d,1)~= size(d,2) || size(d,1) ~= size(d,3)
                disp('input matrix must be square!')
                return
            end
            d_size = size(d,1);
            obj.max_i = d_size-1;
            yzs = [];
            for i = 1:d_size
                yzs = [yzs Spline2D(reshape(d(i,:,:),d_size,d_size))];%append(spline1D.Spline1D(d[i,:]))
            end
            xs = [];
            cx = 0;
            cy = 0;
            while cx <= (obj.max_i+0.01)
                if cx > obj.max_i
                    cx = obj.max_i;
                end
                cy = 0;
                while cy <= (obj.max_i+0.01)
                    if cy > obj.max_i
                        cy = obj.max_i;
                    end
                    xv = zeros(d_size,1);
                    for i = 1:d_size
                        xv(i) = yzs(i).f(cy,cx);
                    end
                    xs = [xs Spline1D(xv)];
                    cy = cy + 1/3;
                end
                cx = cx + 1/3;
            end
            obj.coeff = zeros(obj.max_i, obj.max_i,obj.max_i,64);
            A = zeros(64,64);
            for i = 1:4
                dx = (i-1)/3;
                for j = 1:4
                    dy = (j-1)/3;
                    for k = 1:4
                        dz = (k-1)/3;
                        for l = 1:4
                            for m = 1:4
                                for n = 1:4
                                    A((i-1)*16+(j-1)*4+k,(l-1)*16+(m-1)*4+n) = dx^(l-1)*dy^(m-1)*dz^(n-1);
                                end
                            end
                        end
                    end
                end
            end

            
            b = zeros(64,1);
            row_size = 3*obj.max_i+1;
%             for i = 1:(d_size-1)
%                 for j = 1:(d_size-1)
            for i = 1:obj.max_i
                for j = 1:obj.max_i
                    for k = 1:obj.max_i
                        for m = 1:4
                            for n = 1:4
                                sp = xs((i-1)*3*row_size+(j-1)*3+(m-1)*row_size+n);
                                for o = 1:4
                                    cx = k-1+(o-1)/3;
                                    b((m-1)*16+(n-1)*4+o) = sp.f(cx);
                                end
                            end
                        end
                        x = A\b;
                        obj.coeff(i,j,k,:)=reshape(x(:),[1,1,1,64]);
                    end
                end
            end
            
        end
                   
        

        function yval = f(obj,z,y,x)
%             y=1;
%             x=1;
            [ix, x_diff] = roundAndCheck(x,obj.max_i);
            [iy, y_diff] = roundAndCheck(y,obj.max_i);
            [iz, z_diff] = roundAndCheck(z,obj.max_i);
            if ix == -1|| iy == -1 || iz == -1
                yval = 0;
                return
            end
            yval = 0;
            for i = 1:4
                for j = 1:4
                    for k = 1:4
                        yval = yval+obj.coeff(ix,iy,iz,16*(i-1)+4*(j-1)+k)*x_diff^(i-1)*y_diff^(j-1)*z_diff^(k-1);
                    end
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