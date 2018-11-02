function res = symmetrize(img)
s1 = cat(1,img,img(end:-1:0,:));
res = cat(2,s1,s1(:,end:-1:0));
