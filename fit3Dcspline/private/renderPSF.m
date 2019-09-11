function img=renderPSF(coeff, cor, Npixels)
%par= x,y,z,N,bg

Nfits = size(cor,1);
spline_xsize = size(coeff,1);
spline_ysize = size(coeff,2);
spline_zsize = size(coeff,3);
off = floor(((spline_xsize+1)-Npixels)/2);
data = zeros(Npixels,Npixels,Nfits,'single');

t=tic;
for kk = 1:Nfits
    xcenter = cor(kk,1);
    ycenter = cor(kk,2);
    zcenter = cor(kk,3);
    
    xc = -1*(xcenter - Npixels/2+0.5);
    yc = -1*(ycenter - Npixels/2+0.5);
    zc = zcenter - floor(zcenter);
    
    xstart = floor(xc);
    xc = xc - xstart;
    
    ystart = floor(yc);
    yc = yc - ystart;
    

    zstart = floor(zcenter);
    
   
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2(single(xc),single(yc),single(zc));
    
    for ii = 0:Npixels-1
        for jj = 0:Npixels-1
             temp = fAt3Dj_v2(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);
             model = temp*cor(kk,4)+cor(kk,5);
             data(ii+1,jj+1,kk)=model;
        end
    end
    if toc(t)>1
        disp(kk/Nfits)
        t=tic;
    end
    
end
img=data;
