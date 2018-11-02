function mexSMAP
if isdeployed
    return
end
mexdir='shared/mexfiles/';
currentdir=pwd;
allf=dir([mexdir '*.c*']);
cd(mexdir)
for k=1:length(allf);
    try
        disp(['mexing ' allf(k).name])
        mex([ allf(k).name])
    catch
        disp(['error in ' allf(k).name])
    end
end
cd(currentdir)