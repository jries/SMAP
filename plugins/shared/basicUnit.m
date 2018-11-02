classdef basicUnit < handle

    properties
        %define class properties if needed
        capDepth
        capShape
        bottomDepth
        bottomShape
        radius
        totalDepth
        basicComponent
        basicComponent3D
    end
    methods
        function obj=basicUnit(capDepth, capShape, bottomDepth, bottomShape, radius)
            % basicUnit(25,'ellipse', 80, 'bar', 25)
            obj.capDepth = capDepth;
            obj.capShape = capShape;
            obj.radius = radius;
            obj.bottomDepth = bottomDepth;
            obj.totalDepth = capDepth + bottomDepth;
            obj.bottomShape = bottomShape;
            obj.basicComponent = mkBasicCom(capDepth, capShape, bottomDepth, bottomShape, radius);
            [x,y,z] = cylinder(obj.basicComponent.x, 1000);
            basicComponent3D = table(x(:),y(:),z(:)*obj.totalDepth);
            basicComponent3D.Properties.VariableNames = {'x' 'y' 'z'};
            obj.basicComponent3D = basicComponent3D;
            
        end
        
        function plot(obj, mode)
            switch mode
                case '2d'
                    plot([obj.basicComponent.x -obj.basicComponent.x(end:-1:1)], [obj.basicComponent.y obj.basicComponent.y(end:-1:1)])
                case '3d'
                    scatter3(obj.basicComponent3D.x, obj.basicComponent3D.y, obj.basicComponent3D.z)
                case '3dPx'
                    scatter3(obj.basicComponent3DPx.x, obj.basicComponent3DPx.y, obj.basicComponent3DPx.z)
            end
        end
        
        function time = getTime(obj, v, t)
            time = obj.basicComponent3DPx;
            time.z = obj.basicComponent3DPx.z - obj.totalDepth + v * t;
        end
        
        function fullPx = fullPixel(obj)
            [x,y,z] = meshgrid(round((-obj.radius):obj.radius), round((-obj.radius):obj.radius), round(0:obj.totalDepth));
            fullPx = table(x(:),y(:),z(:));
            fullPx.Properties.VariableNames = {'x' 'y' 'z'};
        end
        
        function ESC = convert2EScoordinates(obj, zmt)
            ESC = obj.basicComponent3D;
            ESC.z = ESC.z+zmt;
        end
    end
    methods(Static)
        function  keept = pointsInCy(X, Y, Z, x, y, z)
            IdxXR = X >= 0;
            IdxXL = X <= 0;
            xRefR = griddata(Y(IdxXR),Z(IdxXR),X(IdxXR), y, z, 'cubic');
            xRefL = griddata(Y(IdxXL),Z(IdxXL),X(IdxXL), y, z, 'cubic');
            IdxYR = Y >= 0;
            IdxYL = Y <= 0;
            yRefR = griddata(X(IdxYR),Z(IdxYR),Y(IdxYR), x, z, 'cubic');
            yRefL = griddata(X(IdxYL),Z(IdxYL),Y(IdxYL), x, z, 'cubic');
            keeptByx = x <= xRefR & x >= xRefL ;
            keeptByy = y <= yRefR & y >= yRefL ;
            keeptByy(y<5&y>-5)=keeptByx(y<5&y>-5);
            keeptByx(x<5&x>-5)=keeptByy(x<5&x>-5);
            keept = keeptByx;
        end
    end
end

function  basicComponent = mkBasicCom(capDepth, capShape, bottomDepth, bottomShape, radius)
    %% make a surface component for making a cylinder
    switch capShape
        case 'ellipse'
        a=capDepth; % horizontal radius
        b=radius; % vertical radius
        xO=0; % ellipse centre coordinates
        yO=bottomDepth;
        t=-pi:0.01:pi;
        x=xO+a*cos(t);
        y=yO+b*sin(t);
        y(y == max(y)) = round(max(y));
        x(y == max(y)) = 0;
        idxKept = x<=xO&y>yO;
        cx = x(idxKept);
        cy = y(idxKept);
    end
   
    if bottomDepth > 0
        switch bottomShape
            case 'bar'
                py = 0:1:bottomDepth;
                px = repelem(-radius, numel(py));
                x = [px cx];
                y = [py cy];
        end
    else
        x = cx;
        y = cy;
    end
    basicComponent = [];
    xy = unique([x;y]','rows','stable');
    yq = 0:(capDepth+bottomDepth);
    xq = interp1(xy(:,2),xy(:,1),yq,'cubic');
    basicComponent.x = xq;
    basicComponent.y = yq;
end

