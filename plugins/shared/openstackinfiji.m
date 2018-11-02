function openstackinfiji(obj,imout3,title)
obj.status('open image in fiji. Large stacks take some time.');
disp('increase the Java heap memory (Matlab -> Preferences -> General) if problems occur duirng opening')
drawnow
    imout=uint8(imout3/max(imout3(:))*(2^8-1));
    s=size(imout3);
    if length(s)>3  
        imout=permute(imout,[1,2,4,3]);
    end
    ijm=openfiji(obj);
    img=copytoImagePlus(imout);
    img.show;
%     ijm.show('imout');
%     mij.createColor(title,imout,true);
    obj.status('done opening stack in fiji.');
end