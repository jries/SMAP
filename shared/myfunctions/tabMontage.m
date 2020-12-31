function tabMontage(nrow)
    f = gcf;
    tabGrp = f.Children;
    nTabs = size(tabGrp.Children,1);

    for k = 1:nTabs
        tabGrp.SelectedTab = tabGrp.Children(k);
        pause(0.5)
        oneFrame = getframe(f);
        if k==1
            frameStack = uint8(zeros([size(oneFrame.cdata) nTabs]));
        end
        frameStack(:,:,:,k) = oneFrame.cdata;
    end
    ncol = ceil(nTabs/nrow);
    figure; montage(frameStack, 'Size', [nrow ncol])
end