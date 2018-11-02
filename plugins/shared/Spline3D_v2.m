classdef Spline3D_v2<handle
    properties
        max_ix;
        max_iy;
        max_iz;
        coeff;
    end
    methods
        function obj = Spline3D_v2(d)
%             if size(d,1)~= size(d,2) 
%                 disp('input matrix must be square in xy!')
%                 return
%             end
            d_sizex = size(d,1);
            d_sizey = size(d,2);
            d_sizez = size(d,3);
            obj.max_ix = d_sizex-1;
            obj.max_iy = d_sizey-1;
            obj.max_iz = d_sizez-1;
%             xys = [];
%             for i = 1:d_sizex
%                 yz = [xys Spline2D(d(:,:,i))];
%             end
%             zs = [];
%             cx =0;
%             cy = 0;
%             while cx<=(obj.max_ix+0.01)
%                 if cx>obj.max_ix
%                     cx = obj.max_ix;
%                 end
%                 cy = 0;
%                 while cy <= (obj.max_iy+0.01)
%                     if cy>obj.max_iy
%                         cy = obj.max_iy;
%                     end
%                     zv = zeros(d_sizez,1);
%                     for i = 1:d_sizez
%                         zv = xys(i).f(cx,cy);
%                     end
%                     zs = [zs Spline1D(zv)];
%                     cy = cy +1/3;
%                 end
%                 cx = cx+1/3;
%             end
%             obj.coeff = zeros(obj.max_)     
            
                 
            for i = d_sizex:-1:1
                s1=Spline2D_v3(reshape(d(i,:,:),d_sizey,d_sizez));
                yzs(i) = s1;%append(spline1D.Spline1D(d[i,:]))
%                 yzs = [yzs Spline2D_v3(reshape(d(i,:,:),d_sizey,d_sizez))];%append(spline1D.Spline1D(d[i,:]))
            end
            maxx=obj.max_iy*obj.max_iz*9;
            xs(maxx) =Spline1D;
            xsind=1;
            %need check cy or cz first
%             cy = 0;
            cz = 0;
            while cz <= (obj.max_iz+0.01)
                if cz > obj.max_iz
                    cz = obj.max_iz;
                end
                cy = 0;
                while cy <= (obj.max_iy+0.01)
                    if cy > obj.max_iy
                        cy = obj.max_iy;
                    end
                    xv = zeros(d_sizex,1);
                    for i = 1:d_sizex
                        xv(i) = yzs(i).f(cy,cz);
                    end
                    s1=Spline1D(xv);
%                     xs = [xs s1];
                    xs(xsind)=s1;
                    xsind=xsind+1;
                    cy = cy + 1/3;
                end
                cz = cz + 1/3;
            end
            xs(xsind:end)=[];
            obj.coeff = zeros(obj.max_ix, obj.max_iy,obj.max_iz,64);
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
            rowy_size = 3*obj.max_iy+1;
            rowz_size = 3*obj.max_iz+1
%             for i = 1:(d_size-1)
%                 for j = 1:(d_size-1)
            for i = 1:obj.max_iz
                for j = 1:obj.max_iy
                    for k = 1:obj.max_ix
                        for m = 1:4
                            for n = 1:4
                                sp = xs((i-1)*3*rowy_size+(j-1)*3+(m-1)*rowy_size+n);
                                for o = 1:4
                                    cx = k-1+(o-1)/3;
                                    b((m-1)*16+(n-1)*4+o) = sp.f(cx);
                                end
                            end
                        end
                        x = A\b;
                        obj.coeff(k,j,i,:)=reshape(x(:),[1,1,1,64]);
                    end
                end
            end
            
        end
                   
        

        function yval = f(obj,z,y,x)
%             y=1;
%             x=1;
            [ix, x_diff] = roundAndCheck(x,obj.max_ix);
            [iy, y_diff] = roundAndCheck(y,obj.max_iy);
            [iz, z_diff] = roundAndCheck(z,obj.max_iz);
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