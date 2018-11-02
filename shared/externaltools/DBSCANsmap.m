function [DBSCANmat, DBSCANtable, eps] = DBSCANsmap(pos, k,Eps,plotaxis,setstatus)

%==========================================================================
%                              FUNCTION
% Performs DBSCAN analysis for cluster identification on molecule position
% file.
%
% Modified version of dbscan.m written by Michal Daszykowski,
% Department of Chemometrics, Institute of Chemistry, 
% The University of Silesia
% http://www.chemometria.us.edu.pl
%
%==========================================================================
%                             INPUT/OUTPUT
% Input: 
%   -pos: .x, .y, .z (only for one channel). or xnm ynm znm
%   -CroppedPos: Position file produced by GSDLoadCrop
%   -Ch1: Which channel in CroppedPos to analyse
%   -k:minimum number of objects in a neighbourhood of an object (minimum
%      cluster size)
%   -Eps: neighbourhood radius, if not known avoid this parameter or put []
%   -plotClusters: generate an image with each cluster circled, 0 = no
%                  image, 1 = geneate image. Unclustered pixels in grey;
%                  clustered are colour-coded
%   -fName: name for saved images
%
% Output: 
%   -DBSCANmat: Structure containing:
%       -.k: k value entered by user
%       -.Eps: epsilon value used for analysis
%       -.Noise: x/y/z array of non-clustered (noise) molecules in the image
%       -.Cluster(ii): structure containing, where ii = cluster number:
%           -.Pos: x/y/z array of x/y/z positions of molecules in the
%                  cluster
%           -.edgePts: for 2D images, provides x/y coordinates of the
%                      bounding edge of the cluster; for 3D images provides
%                      vertice information for use with trisurf
%   -DBSCANtable - n x 5 table:
%       -Columns 1-3: x/y/z molecule positions
%       -Column 4: cluster ID number to which the molecule belongs. -1 =
%                  noise (unclustered) molecule
%       -Column 5: molecule type. core = 1, border = 0, noise/unclustered
%                  = -1)
%   -eps: Epsilon value used (either user-input 'Eps', or determined
%         automatically
%
%==========================================================================
%                             CITATION
%
% This script is provided as a supplemental material in:
%
%   Fabiana A. Caetano, Brennan S. Dirk, Joshua H.K. Tam, P. Craig 
%       Cavanagh, Maria Goiko, Stephen S.G. Ferguson, Stephen H. Pasternak,
%       Jimmy D. Dikeakos, John R. de Bruyn, Bryan Heit. MIiSR: Analysis of 
%		Molecular Interactions in Super-Resolution Imaging Enables the Study
%		of Protein Interactions, Dynamics and Formation of Multi-protein 
%		Structures. 2015 PLoS Computational Biology
% 
% Please reference this paper in any publications which use this script for 
% analysis.
%==========================================================================
%% Check inputs & load file

%Extract x/y/z
fn=fieldnames(pos);
nel=length(pos.(fn{1}));
if isfield(pos,'z') || isfield(pos,'znm')
    x=zeros(nel,3);
else
    x=zeros(nel,2);
end

if isfield(pos,'x')
    x(:,1)=pos.x;
else
    x(:,1)=pos.xnm;
end
if isfield(pos,'y')
    x(:,2)=pos.y;
else
    x(:,2)=pos.ynm;
end
if isfield(pos,'z')
    x(:,3)=pos.z;
elseif isfield(pos,'znm')
    x(:,3)=pos.znm;
end


x2 = x; %preserve original x-variable

[m,n] = size(x); %array size

%Define epsilon, if needed
if nargin < 3 || isempty(Eps)
    Eps=((prod(max(x)-min(x))*k*gamma(.5*n+1))/(m*sqrt(pi.^n))).^(1/n);
end

%% Perform DBSCAN analysis
eps = Eps;
 x=[[1:m]' x];
 xs=sortrows(x,2);
type=zeros(1,m);
no=1;
touched=false(1,m);
countershowo=0;
class=zeros(1,m);
inds1=1;
timerv=tic;
for i=1:m
    if touched(i)==0;
       xph=xs(i,2);
       inds1=find(xs(inds1:end,2)>=xph-eps,1,'first')+inds1-1;
       inds2=find(xs(inds1:end,2)>=xph+eps,1,'first')+inds1-1;
       ob=xs(i,:);
       [indi,nind]=dist2(ob(2:n),xs(:,2:n),inds1,inds2,Eps);
%        ind=find(D<=Eps);
%        indi=D<=Eps;
%        nind=sum(indi);
    
       if nind>1 && nind<k+1       
          type(i)=0;
          class(i)=0;
       end
       if nind==1
          type(i)=-1;
          class(i)=-1;  
          touched(i)=1;
       end

       if nind>=k+1; 
          type(i)=1;
          class(indi)=max(no);
%           class(indi)=ones(nind,1)*max(no);
          while nind>0%~isempty(ind)
%                 ob=x(ind(1),:);
                indi1=find(indi,1,'first');
                if isempty(indi1)
                    nind=0;
                    continue
                end
                ob=xs(indi1,:);
                touched(indi1)=1;
                indi(indi1)=false;
                
%                 inds3=inds2;
%                 inds1b=inds2;
                while(touched(inds1))&&inds1<m
                    inds1=inds1+1;
                end
                while(xs(inds2,2)<ob(2)+eps)&&inds2<m
                    inds2=inds2+1;
                end
%                 inds3=find(xs(inds2:end,2)>=ob(2)+eps,1,'first')+inds2-1;
                [ieps,neps]=dist2(ob(2:n),xs(:,2:n),inds1,inds2,Eps);
%                 i1=find(D<=Eps);
%                 ieps=(D<=Eps);
%                 neps=sum(ieps);
     
                if neps>1
                    
                   class(ieps)=no;
                   if neps>=k+1;
                      type(ob(1))=1;
                   else
                      type(ob(1))=0;
                   end

                   it=ieps&~touched;
                   
                   indi(it)=true;
                   nind=sum(indi);
%                     ind=[ind find(it)];
%                    class(it)=no;
                   touched(it)=true;
%                    for ix=1:length(i1)
%                        if touched(i1(ix))==0
%                           touched(i1(ix))=1;
%                             ind=[ind i1(ix)];   
%                            class(i1(ix))=no;
%                        end                    
%                    end
                   

                   
                end
          end
          no=no+1; 
       end
    end
  
   if toc(timerv)>10
       countershow=sum(touched);
       setstatus(['DBSCAN: clustering localization ' num2str(countershow) ' of ' num2str(m)]);
       drawnow;
%        countershowo=countershow;
       timerv=tic;
   end
end
[~,sortback]=sort(xs(:,1));
class=class(sortback);
i1=find(class==0);
class(i1)=-1;
type(i1)=-1;

%% Generate Output
DBSCANmat.k = k;
DBSCANmat.Eps = eps;
DBSCANmat.Noise = x2(class==-1,:);
if max(class) == -1
    display ('No clusters identified.');
else
    for ii=1:max(class)
        DBSCANmat.Cluster(ii).Pos = x2(class==ii,:);
        try
            if n == 2 %2D
                edgePts = convhull (DBSCANmat.Cluster(ii).Pos);
                DBSCANmat.Cluster(ii).edgePts = DBSCANmat.Cluster(ii).Pos(edgePts,:);
            else %3D
                DBSCANmat.Cluster(ii).edgePts = convhull (DBSCANmat.Cluster(ii).Pos);
            end
        catch
            DBSCANmat.Cluster(ii).edgePts = NaN;
        end
    end
end

if n == 2
    x2(:,3) = 0;
end

DBSCANtable(:,1:3) = x2;
DBSCANtable(:,4) = class;
DBSCANtable(:,5) = type;

if  nargin>3&&~isempty(plotaxis)&&isvalid(plotaxis)
    plotClusters=true;
else
    plotClusters=false;
end
%% Plot image
if plotClusters && max(class) == -1
%     h = figure(1);
    scatter (DBSCANmat.Noise(:,2), DBSCANmat.Noise(:,1), 1, [0.2 0.2 0.2],'Parent',plotaxis);
    axis(plotaxis,'equal'); 
    title(['DBSCAN, no clusters detected. k = ' num2str(k), ', ' char(949) ' = ' num2str(Eps)],'Parent',plotaxis);
elseif plotClusters
    clustercols = lines(max(class));

    scatter (DBSCANmat.Noise(:,2), DBSCANmat.Noise(:,1), 1, [0.7 0.7 0.7],'Parent',plotaxis);
    plotaxis.NextPlot='add';
    for ii=1:max(class)
        scatter(DBSCANmat.Cluster(ii).Pos(:,2), DBSCANmat.Cluster(ii).Pos(:,1), 1, clustercols(ii,:),'Parent',plotaxis);
    end
    axis (plotaxis,'equal');
    title (['DBSCAN, ' num2str(max(class)) ' clusters detected. k = ' num2str(k), ', ' char(949) ' = ' num2str(Eps)],'Parent',plotaxis);
end %end plotting

end


%% Sub-functions
function [indi,N]=dist2(i,x,ind1,ind2,eps)

% function: [D]=dist2(i,x)
%
% Calculates the Euclidean distances between the i-th object and all objects in x	 

% [m,n]=size(x);
% Do=sqrt(sum((((ones(m,1)*i)-x).^2)'));
indi=false(1,length(x));

 r=ind1:ind2;
D=sqrt(((x(r,1)-i(1)).^2+(x(r,2)-i(2)).^2))<=eps;

N=sum(D);
indi(r)=D;



% if n==1
%    D=abs((ones(m,1)*i-x))';
% end

end % end dist2