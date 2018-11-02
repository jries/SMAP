function [imag,info,fileout]=mytiffreadgeneral(file,frames,opt)

%info: gain, size, all...
[pathstr, name, ext]=fileparts(file);


if isempty(name)
    df=myfastdir(pathstr,  'img*.tif');
    whos df
    name=df(1);
    file=[pathstr '/' name];
end
if strfind(name,'img_00') %imageJ micro-manager file
    mm=1;
else
    mm=0
end
    dind=findstr(name,'_');
   
%     ss=[pathstr '/' name(1:dind(1)) '*' name(dind(2):end)  ext]
% ss=[pathstr '/' name(1:dind(1)) '*' ext] %load all files in directory
%     allfiles=dir(ss)
searchx=[name(1:dind(1)) '*' ext];
allfiles=myfastdir(pathstr,searchx);
    fileout.names=allfiles;
    if frames==0
        frames=1:length(allfiles);
    end
    fileout.pfad=[pathstr '/'];
    fileout.frames=frames;
    fileout.numframes=length(allfiles);
    

    finfo=imfinfo(file);
    [p n x]=fileparts(file);
    
    
    metafile=([ p filesep 'metadata.txt']);
    fid=fopen(metafile,'r');
    if fid>0 %MM metadata present
        fclose(fid);
        minfo=fileread(metafile);
        info=minfoparse(minfo);

%        info.actualconv
    %further process info
        if ~isnan(info.actualconv)   
            info.conversion=info.actualconv;
        else
            if strfind(minfo,'Evolve')
            if ~isempty([strfind(info.port,'EM') strfind(info.port,'Multiplication')])
                disp('Evolve EM')
                switch info.gain %for QuantEM
                    case 1
                        info.conversion=3.9; %preliminary..
                    case 2
                        info.conversion=2.3; %preliminary.. 
                    case 3
                        info.conversion=3; %preliminary.. 
                end
                info.emgainnom=info.emgain;
        %         x=log(info.emgainnom);
        %         y = 0.01497*x^3 - 0.209*x^2 + 1.528*x - 0.683;
        %         y = 0.011859*x^3 - 0.15177*x^2 + 1.2036*x - 0.040978;
        %         
        %         info.emgain=exp(y);
            else
                    switch info.gain %for QuantEM
                    case 1
                        info.conversion=4; %preliminary..
                    case 2
                        info.conversion=2; %preliminary.. 
                    case 3
                        info.conversion=1; %preliminary.. 
                    end
                    info.emgain=1;
                    info.offset=500;
            end
            info.cam_pixelsize_um=0.1;%38;  %depends on settings, camera
            elseif strfind(minfo,'Andor sCMOS')
                %only EM-Gain, 10 MHz, baseline clamp on
                 info.offset=100;

                 info.conversion=.666; %preliminary..
                 info.actualconversion=.666;
                 info.cam_pixelsize_um=0.096;  %2x2 binning, f=400
                  info.port='CCD';
                  disp('Andor  sCMOS')
            elseif strfind(minfo,'iXon Ultra')
                info.offset=200;
                switch info.port
                    case 'Conventional'
                        info.port='Conventional';
                        info.emgain=1;
                    switch info.readoutrate
                        case '3.000 MHz'
                            mode=5;
                        case '1.000 MHz'
                            mode=6;
                        case '0.080 MHz'
                            mode=7;
                    end
                        
                    case 'Electron Multiplying'
                        info.port='EM';
                    switch info.readoutrate
                        case '17.000 MHz'
                            mode=1;
                        case '10.000 MHz'
                            mode=2;
                        case '5.000 MHz'
                            mode=3;
                        case '1.000 MHz'
                            mode=4;
                    end
                end
                gainmatrix=[14.61 9.3 4.81
                    15.46 7.7 4.2
                    16.73 7.87 4.2
                    16.77 8.04 4.1
                    3.83 3 1.4
                    3.87 2.95 1.37
                    3.84 2.96 1.36];
                info.conversion=gainmatrix(mode,info.gain);
                info.cam_pixelsize_um=0.138;

             elseif strfind(minfo,'Andor')
                %only EM-Gain, 10 MHz, baseline clamp on
                info.offset=100;
                if strfind(info.port,'Off')
                    info.emgain=1;
                    switch info.gain
                        case 5.1
                            info.conversion=10.39; %%%Markus
                        otherwise
                            info.conversion=0;
                    end
                else
                 
                switch info.gain %Andor from Malte
                    case 4.6
                        info.conversion=12; %preliminary..
                    case 2.3
                        info.conversion=26.3; %preliminary.. 
                    case 1
                        info.conversion=62.2; %preliminary.. 
                    otherwise
                        info.conversion=62.2;
                end
                end
                info.cam_pixelsize_um=0.138;  %depends on settings, camera
            else %default
                disp('dummy 1')
                info.emgain=18;  
                
                info.conversion=1;
                info.offset=4000;
                info.cam_pixelsize_um=0.13;  %depends on settings, camera
            end
       
        

        end
%         if ~strfind(info.readoutrate,'10MHz')
%             info.conversion=info.conversion/3;  %%No EM: careful: 5 MHz. 
%         end
    else %no metafile present: now: assign dummy values. later read sif etc.
        mf=dir([pathstr '/' name(1:dind(1)-1) '*.sif']);
        mff=[pathstr '/' mf.name];
        fid=fopen(mff);
        if fid>0 %Andor SIF file
            fclose(fid);
            info=mysifinfo2(mff);
        else
        disp('here we changed values: Anna')
%             info.conversion=12;
%             info.actualconv=10.9;
%             info.conversion=10.88;
%             info.actualconv=10.88;
            info.actualconv=10;
            info.conversion=10;

            info.readoutrate=1;
            info.exposure=1;
            info.gain=1;
            info.offset=0;
            info.emgain=100;  
            info.chip='not determined';
            info.timediff=1;
            info.port='EM';
            info.roi=[];
            info.cam_pixelsize_um=0.1;  %depends on settings, camera
        end

    end
    


info.Width= finfo(1).Width;
info.Height=finfo(1).Height;


w=finfo(1).Width;
h=finfo(1).Height;
% finfo(1).BitDepth
if finfo(1).BitDepth==16
    vtype='uint16';
elseif finfo(1).BitDepth==8
    vtype='uint8';
end

if nargin==3
    if strcmp(opt,'numf')
        dx=ceil(fileout.numframes/frames);
        frames=ceil(dx/2):dx:fileout.numframes;
    end
end

nframes=length(frames);
entries=h*w*nframes;

% [user,sys] = memory;
% free=sys.PhysicalMemory.Available
free=3e9;
maxmem=free*0.7;
if entries*4>maxmem %4GB max
    nframes=round(maxmem/h/w/4)
    errordlg(['too much memory required, loading only first ' num2str(frames) ' frames'])

end

imag=zeros(h,w,nframes,vtype); 

fmt_s = imformats('tif');

for k=1:nframes
    if ~abs((round(k/10)-k/10))
%         disp(k)
    end
    if frames(k)<=length(allfiles)
    fr=myimread([pathstr '/' allfiles{frames(k)}],fmt_s);
%         imag(:,:,k)=single(fr);
        imag(:,:,k)=(fr);
        
%         figure(1)
%         imagesc(fr)
%         asdf
    else
        break;
    end
end

info.frames=frames;
fileout.frames=frames;
info.size=[h w frames];

% if ~mm
%     info.offset=min(imag(:));
% end
