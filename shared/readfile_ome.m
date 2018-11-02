function img=readfile_ome(file)
try
imgall=bfopen(file);
catch
    disp(['could not open file ' file])
    img=[];
    return
end
sim=size(imgall{1}{1});
numframes=length(imgall{1});
img=zeros(sim(1),sim(2),numframes);
for k=1:numframes
    img(:,:,k)=imgall{1}{k};
end
end