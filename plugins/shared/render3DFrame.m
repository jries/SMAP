function oneFrame=render3DFrame(locs,f)

      p.sgauss=[1.5, 0.8];%smoothing of volume
            p.cutoff=4; % for isosurface, at max(V(:))/cutoff
            p.cutoffdc=[2.2 3]; % for isosurface, at max(V(:))/cutoff
            p.pxSize = 5; % pixel size
            
            mx = [0 200];
            my = [-200 200];
            mz = [-200 200];
            
            x = locs.xnm-1000;
            y = locs.ynm;
            z = locs.znm;
            x(x<0) = -x(x<0);
            
            %layer 1
            inl1=locs.layer==1;
            [v, xb, yb,zb] = myhist3(x(inl1), y(inl1), z(inl1), p.pxSize, mx,my,mz);
            vs=smooth3(v,'gaussian',[1 1 1]*(round(2*p.sgauss(1)/2)*2+1),p.sgauss(1));
              co=max(vs(:))/p.cutoffdc(1);    
%               co=4;
            %layer 2
            
            inl2=locs.layer==2;
            [v2, xb, yb,zb] = myhist3(x(inl2), y(inl2), z(inl2), p.pxSize, mx,my,mz);
            vs2=smooth3(v2,'gaussian',[1 1 1]*(round(2*p.sgauss(2)/2)*2+1),p.sgauss(2));
             codc2=max(vs2(:))/p.cutoffdc(2);  
             codc2=0.25;
            
            clf(f)
            ax = axes(f);
            ax.Color=[0 0 0];
            fc = [1,.75,.65];
%             fc = 'red';
            hiso = patch(ax, isosurface(vs,co),'EdgeColor','none','FaceColor',fc);
            isonormals(vs,hiso);
            hold on
            hcap=patch(ax, isocaps(vs,co),'FaceColor','interp','EdgeColor','none');
            
            
            hiso2=patch(isosurface(vs2,codc2),'EdgeColor','none','FaceColor',[0,1,1]);
           hold on
            isonormals(vs2,hiso2);
            hcap2=patch(ax,isocaps(vs2,codc2),'FaceColor',[0,.7,1],'EdgeColor','none');
%         hcap2=patch(ax,isocaps(vs2,codc2),'FaceColor','interp','EdgeColor','none');

            
            
            axis equal
            xlim([0 range(my)/p.pxSize])
            ylim([0 range(mx)/p.pxSize])
            zlim([0 range(mz)/p.pxSize])
            view(-35,30)
            colormap hot
%             xlim([0,size(vs,1)])            
            lightangle(45,30);
            lighting gouraud
            colormap hot;
            drawnow
            oneFrame = getframe(f);
%             oneFrame = oneFrame.cdata;