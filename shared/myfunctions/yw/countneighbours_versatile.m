function neighboursback=countneighbours_versatile(ref,target,dxyz,region)
    sortmallr=horzcat(ref(:,1),(1:size(ref,1))');
    [sortedmr,sortindr]=sortrows(sortmallr);
    refsorted=ref(sortindr,:);

    sortmallt=horzcat(target(:,1),(1:size(target,1))');
    [sortedmt,sortindt]=sortrows(sortmallt);
    targetsorted=target(sortindt,:);

    if size(ref,2)==3 %3D
        if region==1 %Gauss
            countf=@countneighbours3DGauss2;
        else
            countf=@countneighbours3Dcirc2;
        end
        neighbours=countf(double(refsorted(:,1)),double(refsorted(:,2)),double(refsorted(:,3)),...
                double(targetsorted(:,1)),double(targetsorted(:,2)),double(targetsorted(:,3)),double(dxyz(1)),double(dxyz(end)));   
    else
        if region==1 %Gauss
            countf=@countneighbours2DGauss2;
        else
            countf=@countneighbours2Dcirc2;
        end
        neighbours=countf(double(refsorted(:,1)),double(refsorted(:,2)),...
                double(targetsorted(:,1)),double(targetsorted(:,2)),double(dxyz(1)));
    end
    [~,sortbackind]=sort(sortedmr(:,2));
    neighboursback=neighbours(sortbackind);
end