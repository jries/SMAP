function out=chooseTifImage(filename,P)
%il=imageloaderOME(filename,[],P);
il=imageloaderAll(filename,[],P);
% il=imageloaderTifSimple(filename);
numberOfFrames=il.metadata.numberOfFrames;
allframes=il.getmanyimages(1:numberOfFrames,'mat');
if numberOfFrames==1
    out=makeout(allframes,il);

    return;
end

thisframe=1;

f=figure;
ax=axes(f,'Position',[0.02,0.15,0.65,0.83]);
imagesc(ax,allframes(:,:,thisframe))

fs=get(0,'defaultUicontrolFontSize')*1.5;
rangenum=1:numberOfFrames;

slider=uicontrol('style','slider','Units','normalized','Position',[0.02,0.02,0.5,.05],'Callback',@slidercallback,'SliderStep',[1 10]/numberOfFrames,'Max',numberOfFrames,'Min',1,'Value',1);
number=uicontrol('style','edit','Units','normalized','String',1,'Position',[0.52,0.02,0.13,.05],'Callback',@numbercallback,'FontSize',fs);
mode=uicontrol('style','popupmenu','Units','normalized','String',{'current','all in range','mean','MIP'},'Position',[0.7,0.8,0.28,0.05],'FontSize',fs);
rtext=uicontrol('style','text','Units','normalized','String','range:','Position',[0.7,0.7,0.28,0.05],'FontSize',fs);
range=uicontrol('style','edit','Units','normalized','String',['1:' num2str(numberOfFrames)],'Position',[0.7,0.6,0.28,0.05],'FontSize',fs);

ok=uicontrol('style','pushbutton','Units','normalized','String','ok','Position',[0.7,0.02,0.28,0.1],'FontSize',fs,'Callback',@ok_callback);

ok=false;
waitfor(f)

switch modeselection
    case 'current'
        outim=allframes(:,:,thisframe);
    case 'all in range'
        outim=allframes(:,:,rangenum);
    case 'mean'
        outim=mean(allframes(:,:,rangenum),3);
    case 'MIP'
        outim=max(allframes(:,:,rangenum),[],3);
end
out=makeout(outim,il);



    function slidercallback(a,b)
        number.String=(num2str(round(a.Value)));
        thisframe=round(a.Value);
        imagesc(ax,allframes(:,:,thisframe));
    end

    function numbercallback(a,b)
        thisframe=round(str2double(a.String));
        slider.Value=thisframe;
        imagesc(ax,allframes(:,:,thisframe));
    end
    function ok_callback(a,b)
        modeselection=mode.String{mode.Value};
        rangenum=str2num(range.String);
        ok=true;
        close(f);
    end

    function imout=makeout(image,il)
        sim=size(image);  
        if length(sim)<3
            sim(3)=1;
        end
        for k=sim(3):-1:1
            imout(k).image=image(:,:,k); 
            imout(k).info.Width=sim(1);
            imout(k).info.Height=sim(2);
            imout(k).info.roi=il.metadata.roi;
            imout(k).info.name=strrep(il.file,'.tif',['_' num2str(k) '.tif']);
        end
    end

end



