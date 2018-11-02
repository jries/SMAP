classdef AnalyzeSPT<interfaces.DialogProcessor
    %
    methods
        function obj=AnalyzeSPT(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson','cam_pixelsize_nm'};
            obj.showresults=true;
        end
        function out=run(obj,p)
            
            out=analyzei(obj,p);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end
function out=analyzei(obj,p)
out=[];
[locs,indin]=obj.locData.getloc({'xnm','ynm','znm','frame','track_id','diffusionCoefficient'},'layer',1,'position','roi','grouping','ungrouped');
intrack=find(locs.track_id>0);
trackid=locs.track_id(intrack);
%maybe filter with minimum length here: histogram 1:N
hc=histcounts(trackid,1:max(trackid)+1);
minlen=10;
goodids=(hc>minlen);
longtrack=goodids(trackid);
inindin=intrack(longtrack);
loct.x=locs.xnm(intrack(longtrack));
loct.y=locs.ynm(intrack(longtrack));
loct.frame=locs.frame(intrack(longtrack));
loct.track_id=locs.track_id(intrack(longtrack));

[~,sortid]=sort(loct.track_id);

loct.x=loct.x(sortid);
loct.y=loct.y(sortid);
loct.frame=loct.frame(sortid);
loct.track_id=loct.track_id(sortid);

if ~isempty(locs.znm)
    loct.z=locs.znm(intrack(longtrack));
    loct.z=loct.z(sortid);
end


p.timediff=30;%XXX
p.mintracklength=minlen;
img=obj.getPar('sr_image');
ax=obj.initaxis('tracks overlay');
imagesc(ax,img.rangex*1000,img.rangey*1000,img.image)
hold on
statistics=plottracks(loct,p,2);
hold off

switch p.analysismode.selection
    case 'statistics'
        getstatistics(statistics,p,obj)

    case 'interactive track exploration'
        avim=sum(img.image,3);
        c = max( avim(:) );
        avim = avim / c;
        avim2 = zeros([size(avim) 3]);
        avim2(:,:,1) = avim;
        avim2(:,:,2) = avim;
        avim2(:,:,3) = avim;
        plotTrackDiffusion(loct,p,img.rangex*1000,img.rangey*1000,avim2,-2,2);
    case 'diffusion maps'
        makediffusionmaps(statistics,loct,p,obj)
    case 'grid based diffusion coefficients'
        	f = figure(34)
% 			set(gcf, 'Position', get(0,'Screensize')); % maximize
			
			% get minD (GUI)
			minD = 10^p.cutoffimmobile;
            statt=statistics;
			% consider only tracks with log10(diffusion) >= minD
			goodt = (statt.tracks(:,8) > 0) & (log10(statt.tracks(:,8)) >= minD);
			t = statt.tracks(goodt,:);
			
			% Listbox
			u1 = uicontrol('Style', 'popup',...
				'String', 'Track Average|weigthed OLS|Robust Fit',...
				'Position', [10 800 150 100],...
				'Callback', {@plotGridCells, t, par, pos, f});
			
			% plot grid cells
			usedD = plotGridCells(0,0,t,par,pos,f);
			
			% plot histogram of grid cell diffusion
			subplot(2,2,2)
			hist(usedD)
			title('Diffusion Histogram (all cells, track-based)')
end

%write diffusion coefficient to localization data.
if p.saveD
    D=statistics.tracks/4; %slope=4D
    if 1%isempty(locs.diffusionCoefficient)
        idall=zeros(size(indin),'single')-1;
    else
        idall=obj.locData.loc.diffusionCoefficient;
    end
    findin=find(indin);
    idall(findin(inindin))=D;
    obj.locData.setloc('diffusionCoefficient',idall);
    out=[];
    obj.locData.regroup;
    obj.setPar('locFields',fieldnames(obj.locData.loc))
end


end
        
function getstatistics(statt,p,obj)
% ax=obj.initaxis('tracks overlay');
% img=obj.getPar('sr_image');
% imagesc(ax,img.rangex*1000,img.rangey*1000,img.image)
% hold on
% % statt=plottracks(loct,p,2);
% hold off

    cutoffimmobile=p.cutoffimmobile;		
%     statt = getD(tracks,par,0);
    slo = statt.slo;
    off = statt.off;
    Ds = slo/4;

ax=obj.initaxis('D hist');
maxD=myquantile(Ds,.995);
minD=myquantile(Ds,.005);
step=10^floor(log10(maxD)-1);
    n = minD:step:maxD;
    h = hist(Ds,n);
    bar(n,h);
    parmhat = lognfit( Ds(Ds>0) );
    hold on
    pp=lognpdf(n,parmhat(1),parmhat(2));
    pp=pp/sum(pp)*sum(h(n>=0) );
    plot(n,pp,'r')
    hold off

    D = median(slo)/4; %2D
    title(['Dmed=' num2str(D) ' , logn mu=' num2str(    exp( parmhat(1) ))])
    xlabel('D µm^2/s')

    s = sqrt(median(off));
ax=obj.initaxis('offset');
    hist((off),50)
    title(['sqrt(median(off))=' num2str(s)])
    xlabel('offset^2 (nm^2)')

    ax=obj.initaxis('logD');
    indb = slo<=0;
    Di = slo/4;
    Di(indb) = 1e-8;
    Dlg = log10(Di);

    nlg = log10(max(1e-8,minD)):0.1:log10(maxD);
    xlabel('log10 D')

    indnorm = find(nlg>cutoffimmobile,1,'first');

    h = hist(Dlg,nlg);
%     alldat(ind).nDl = nlg;alldat(ind).hDl=h;
    bar(nlg,h);
    hold off
ax=obj.initaxis('cumulateive sum');
    Dlgs = cumsum(h);
    Dlgs = Dlgs-Dlgs(indnorm);
    Dlgs = Dlgs/Dlgs(end);
    plot(nlg,Dlgs);
    xlabel('cumulative log D')
%     alldat(ind).nDc = nlg;alldat(ind).hDc=Dlgs;
ax=obj.initaxis('tracklength');
    nlen = 1:50;
    hlen = hist(statt.lent,nlen);
    bar(nlen,hlen);
    xlabel('tracklength')
    title(['median tracklength = ' num2str(median(statt.lent)) ' total nr of tracks ' num2str(sum(hlen))]);
    out=[];
end

function makediffusionmaps(statt,locs,p,obj)
% minD=myquantile(Ds,.005);
minD = 10^p.cutoffimmobile;
trD=statt.tracks;
goodt=find(trD>10^p.cutoffimmobile);
[x,sortind]=sort(locs.x(goodt));
y=locs.y(goodt(sortind));
D=trD(goodt(sortind));
% 			statt = getD(tracks,par,0);
% 			goodt = statt.tracks(:,5)>1e-4; % only D>0
% 			tr = statt.tracks(goodt,:);
% 			[x,indx] = sort(tr(:,1));
% 			y = tr(indx,2);
% 			D = tr(indx,5);
			
			%define grid
			dg = p.cam_pixelsize_nm; % camera pixels
			xmin = min(x);
			xmax = max(x);
			ymin = min(y);
			ymax = max(y);
			xg = xmin:dg:xmax;
			yg = ymin:dg:ymax;
			
			Dmap = zeros(length(xg)-1,length(yg)-1);
			
			ind1 = 1;
			for ixh=1:length(xg)-1;
				while x(ind1)<xg(ixh)
					ind1=ind1+1;
				end
				
				ind2=ind1;
				
				while x(ind2)<xg(ixh+1)
					ind2=ind2+1;
				end
				
				indh=ind1:ind2-1;
				xposs=x(indh);yposs=y(indh);Dposs=D(indh);
				
				for iyh=1:length(yg)-1;
					ygood=yposs>yg(iyh)&yposs<yg(iyh+1);
					Dg=Dposs(ygood);
					Dm=median(Dg(Dg>minD));
					Dmap(ixh,iyh)=Dm;
				end
            end
			
            image=obj.getPar('sr_image');
            xl=image.rangex;yl=image.rangey;
            avim=sum(image.image,3);
% 			xl=get(par.handles.axes_rec,'XLim')
% 			yl=get(par.handles.axes_rec,'YLim')
% 			
% 			bgim=pos.palmreconstruct.palmreconstruct;
% 			sb=size(bgim);
% 			
% % 			avim=rescaleimage(bgim,par.gamma,par.quantile);
% 			c=max(avim(:));
% 			avim(avim>c)=c;avim=avim/c;
			%     avim2=zeros([size(avim) 3]);
			%     avim2(:,:,1)=avim;avim2(:,:,2)=avim;avim2(:,:,3)=avim;
			%
			%
            Dmap=Dmap';
			minDl=p.cutoffimmobile;
			maxD=log10(myquantile(D,.995));
			
			mask=1-2*avim;
			mask(mask<0)=0;
			figure(25)
			ind0=isnan(Dmap);
			ind0=ind0|isinf(Dmap);
			Dmap(ind0)=0;
			Dmap(Dmap<=0)=1e-8;
			Dl0=log10(Dmap);
			%     Dl(isnan(Dl))=0;
			Dl=Dl0-minDl; %from -3 : 1;
			Dl=Dl/(maxD-minDl);
			Dl(Dl>1)=1;
			
			%     Dl=Dmap/max(Dmap(:));
			%     ind0=Dl<=0;
			ind0=ind0|isnan(Dl);
			Dl(Dl<0)=0;
			Dldraw=Dl0; Dldraw(ind0)=minDl; Dldraw(Dldraw>maxD)=maxD; Dldraw(Dldraw<minDl)=minDl;
			Dldraw(1,1)=maxD;
			ax=obj.initaxis('median diffusion map');
			imagesc(xl,yl,Dldraw)
			colormap jet
			colorbar
			savim=size(avim);
            Dlhr=imresize(Dl,savim([1 2]),'nearest');
            
% 			[x,y]=ind2sub(size(Dl),find(ind0));
			
			zv=0.4
			
			cm=jet(256);
			Drgb=ind2rgb(uint8(Dlhr*255),cm);
			
% 			i1=sub2ind(size(Drgb),x,y,1*ones(length(x),1));
% 			i2=sub2ind(size(Drgb),x,y,2*ones(length(x),1));
% 			i3=sub2ind(size(Drgb),x,y,3*ones(length(x),1));
% 			Drgb(i1)=zv;Drgb(i2)=zv;Drgb(i3)=zv;
% 			
% 			bf=(1000*par.pixelsize)/par.pixrec*dg;
% 			
% 			Db=imresize(Drgb,bf,'nearest');
% 			sDb=size(Db);
% 			
% 			cr=par.roicoords;
% 			ys=round((cr(1)-xl(1))/par.pixrec*1000)
% 			xs=round((cr(2)-yl(1))/par.pixrec*1000)
% 			Dim=zeros(sb(1),sb(2),3)+zv;
% 			Dim(xs+1:sDb(1)+xs,ys+1:sDb(2)+ys,:)=Db;
			
			ax=obj.initaxis('diffusion map in SR');
			
% 			size(Db)
			h=imshow(Drgb);
			
% 			     colormap jet
% 			     colorbar
			hold on
			k=   imshow(0*sum(Drgb,3));
			hold off
			set(k,'AlphaData',mask)
end

function out=plottracks(locs,p,show)

	pfound = max(locs.track_id);
	cols = jet(pfound); 		% default colormap
	off = zeros(pfound,1);		% double vector (all zeros)
	slo = off;
	lent = off;
	indg = false(size(locs.track_id));			% boolean vector (all false)
	ind1 = 1;
	st = length(locs.track_id);
	dt = p.timediff/1000;
	
    p.minlenplot=10;
	to=zeros(st,1);			% add a vector with zeros

	for k = 1:pfound
		
		% increment ind1 as long as
		%  - it is greater then the number of rows in tracks
		%  - the corrisponding row has track-id smaller k
		while( ind1 < st(1) && locs.track_id(ind1) < k )
			ind1 = ind1 + 1;
		end
		
		ind2 = ind1;
		
		% increment ind2 as long as
		%  - it is greater then the number of rows in tracks
		%  - the corrisponding row has track-id equal k
		while( ind2 < st(1) && locs.track_id(ind2) == k)
			ind2=ind2+1;
		end
		
		% take all rows between ind1 and ind2 (-> index always k)
		ind = ind1:ind2 - 1;
		
		% par.analv1 == min L (GUI)
		if ( length(ind) >= p.minlenplot)
			indg(k) = true; 					% mark row as checked
			x = locs.x(ind) ;
			y = locs.y(ind);
			d = zeros( length(x)-1 ,1);			% vector with length of x
			
			% track step size considered for diffusion
			mint = 2;
			maxt = 4;
			
			% d contains the sum of mean squared distances for dt = 2 : 4 (time step size)
			for l = mint:maxt
				d(l) = mean( (x(1:end-l)-x(l+1:end)).^2 + (y(1:end-l)-y(l+1:end)).^2 );
			end
			
			t = dt * (mint:maxt)';
			
			% X(:,1) = 1, X(:,2) = t
			X = [ones(size(t)) t];
			
			% X  *  coeffs  =  msd
			coeffs = X \ d(mint:maxt);
			slo(k) = coeffs(2)/1e6;         % slope  (k)
			off(k) = coeffs(1);         % offset (d)
			lent(k) = length(ind);		% track length
			
			to(ind,end) = coeffs(2);	% add slope to tracks
			
			if (show == 2)	%plot tracks
				plot(x,y,'Color',cols(k,:))
			end
			
			if (show == 1)						%plot msd vs time
				for l = 1:length(x)-1
					d(l) = mean( (x(1:end-l)-x(l+1:end)).^2 + (y(1:end-l)-y(l+1:end)).^2 );
				end
				plot(d)
				hold on
			end
		end
	end

	out.slo = slo(indg);
	out.off = off(indg);
	out.lent = lent(indg);
	out.tracks = to/1e6;

end

% plotTrackDiffusion
%
% Input
% - tracks   		 : [x y time id ...]
% - par	   	 		 : parameters (global variable)
% - imgx,imgy,avim2	 : background image parameters
% - minD,maxD        : min and max diffusion coefficient
function plotTrackDiffusion(locs,par,imgx,imgy,avim2,minD,maxD)

pfound = max(locs.track_id);

	st = length(locs.track_id);
	dt = par.timediff/1000;
    
    colors=jet(255);
	ind1 = 1;

	figure(33)
	clf
	% plot background image
	image(imgx,imgy,avim2);
	
	%Min
	uicontrol('Style', 'slider',...
		'Min',-8,'Max',2,'Value',minD,...
		'Position', [80 50 500 20],...
		'Callback', {@plotTrackDiffusionHelper,locs,par,imgx,imgy,avim2,minD,maxD,1});
	%Max	
	uicontrol('Style', 'slider',...
		'Min',-4,'Max',6,'Value',maxD,...
		'Position', [80 20 500 20],...
		'Callback', {@plotTrackDiffusionHelper,locs,par,imgx,imgy,avim2,minD,maxD,0});
	uicontrol('Style', 'pushbutton', 'String', 'Adjust Colors',...
		'Position', [600 20 100 20],...
		'Callback', 'colormapeditor');
	text = uicontrol('Style','text',...
		'Position',[600 45 100 20],...
		'String','D');

	hold on
	Min = inf;
	Max = -inf;
	usedD = [];
	for k = 1:pfound
		
		% increment ind1 as long as
		%  - it is greater then the number of rows in tracks
		%  - the corrisponding row has track-id smaller k
		while( ind1 < st(1) && locs.track_id(ind1) < k )
			ind1 = ind1 + 1;
		end
		
		ind2 = ind1 + 1;
		
		% increment ind2 as long as
		%  - it is greater then the number of rows in tracks
		%  - the corrisponding row has track-id equal k
		while( ind2 < st(1) && locs.track_id(ind2) == k)
			ind2=ind2+1;
		end
		
		% take all rows between ind1 and ind2 (-> index always k)
		ind = ind1:ind2 - 1;
		
		% par.analv1 == min L (GUI)
		if ( length(ind) >= par.mintracklength) 
			indg(k) = true; 					% mark row as checked
			x = locs.x(ind);
			y = locs.y(ind) ;
			d = zeros( length(x)-1 ,1);			% vector with length of x
			
			% track step size considered for diffusion
			mint = 2;
			maxt = 4;
			
			% d contains the sum of mean squared distances for dt = 2 : 4 (time step size)
			for (l = mint:maxt)
				d(l) = mean( (x(1:end-l)-x(l+1:end)).^2 + (y(1:end-l)-y(l+1:end)).^2 )/1e6;
			end
			
			t = dt * (mint:maxt)';
			
			% X(:,1) = 1, X(:,2) = t
			X = [ones(size(t)) t];
			
			% X  *  coeffs  =  msd
			coeffs = X \ d(mint:maxt);
            coeffs(imag(coeffs)~=0)=0;
% log10(coeffs(2) / 4)
			colormap jet;
			if coeffs(2)>0&&( log10(coeffs(2) / 4) < maxD && log10(coeffs(2) / 4) > minD )
                usedD = [usedD coeffs(2)/4];
				Min = min(coeffs(2)/4,Min);
				Max = max(coeffs(2)/4,Max);
                colcoeff=coeffs(2)/4;colcoeff(colcoeff<0)=0;
                colcoeff=(log10(colcoeff)-minD)/(maxD-minD);
                c=round(colcoeff*255); c(c<1)=1;c(c>255)=255;
                l = plot(x,y,'Color',colors(c,:));
% 				l = color_line(y,x,repmat(coeffs(2)/4,1,length(x)));
				set(l,'ButtonDownFcn',{@plotTrackD coeffs(2)/4 text})
			end
		end
    end
    
	s = (maxD - minD) / 10;
    cr=0:0.1:1;
	h = colorbar;
	set(h,'YTick',cr)
    
% 	set(h,'YLim',[minD maxD])
	set(h,'YTickLabel',{(minD:s:maxD)})
% % 	set(h,'YTickMode','auto')
% 	set(gcf, 'Renderer', 'opengl')

	hold off
	axis equal
	axis([ min(locs.x) max(locs.x) min(locs.y) max(locs.y) ])
	title(['D interval = [' num2str(minD,2) ',' num2str(maxD,2) '] - Median D = ' num2str(median(usedD),3) ' µm^2/s'])
end	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotTrackDiffusionHelper
% updates min and max value and executes plotTrackDiffusion
%
% Input
% - hObj	 	     : slider handle
% - event	 		 : event data
% - tracks   		 : [x y time id ...]
% - par	   	 		 : parameters (global variable)
% - imgx,imgy,avim2	 : background image parameters
% - minVal, maxVal   : time lag interval for MSD
% - newMin           : boolean (1 == new value is new min, 0 == new value is new max)
function plotTrackDiffusionHelper(hObj,event,tracks,par,imgx,imgy,avim2,minD,maxD,newMin)
	if hObj == 0
		plotTrackDiffusion(tracks,par,imgx,imgy,avim2,minD,maxD);
	else 
	    val = get(hObj,'Value');
	    if (newMin)
	        minD = val;
	    else
	        maxD = val;
	    end
	end

    plotTrackDiffusion(tracks,par,imgx,imgy,avim2,minD,maxD)
end
% plotTrackD

% Input
% - hObj : line handle
% - event: event data
% - d	 : diffusion coefficient
% - hTxt : text handle
function plotTrackD(hObj,event,d,hTxt)
	for h = findobj('LineWidth',2)
		set(h,'LineWidth',.5);
	end
	set(hObj,'LineWidth',2);
	set(hTxt,'String',strcat('D = ', num2str(d)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotGridCells:
% Plot grid-based diffusion map
%
% Input
% - hObj.Value: * 1 => track-based average grid cell diffusion (***default***)
%   			* 2 => ensemble-averaged grid cell diffusion (weighted OLS)
%   			* 3 => ensemble-averaged grid cell diffusion (robust OLS)
% - f		  : figure handle
%
% Output
% - usedD	  : Diffusion coefficients for histogram
function usedD = plotGridCells(hObj, event, tracks, par, pos, f)
	if hObj == 0
		mode = 1;
	else
		mode = get(hObj,'Value');
	end

	% x and y = coordinates of tracks for background image
	% x = y axis?, y = x axis?
	c = par.cpix;
	x = ( (c(3):c(4)+1/par.refinement+.5) + 1/par.refinement ) * par.pixelsize;
	y = ( (c(1):c(2)+1/par.refinement+.5) + 1/par.refinement ) * par.pixelsize;

	avim = pos.palmreconstruct.palmreconstruct;
	c = max( avim(:) );
	avim = avim / c;
	avim2 = zeros([size(avim) 3]);
	avim2(:,:,1) = avim;
	avim2(:,:,2) = avim;
	avim2(:,:,3) = avim;
	subplot(2,2,1)

	% plot background image
	image(y,x,avim2);
	hold on

	% plot tracks
	getD(tracks,par,2);

	% grid cell parameters
	dx = par.analv2;
	dy = par.analv2;

	tracks(:,1:2) = tracks(:,1:2) .* par.pixelsize;
	% get grid coordinates
	grid = getGridIndex(tracks,dx,dy);
	gridD = [];
	usedD = [];

	switch mode
		case 1
			% weighted avg (track-based)
			set(f,'Name', 'Track-based diffusion')
			gridD = getAvgGridCellDiffusion(grid,par,dx,dy);
		case 2
			% weighted least squares fit
			set(f,'Name', 'ensemble-averaged diffusion')
			gridD = getEnsembleLSFDiffusion(grid,par,dx,dy,0);
		case 3
			% robust fit
			set(f,'Name', 'robust linear regression')
			gridD = getEnsembleLSFDiffusion(grid,par,dx,dy,1);
    end
    if mode == 1
        title(strcat('D = ',num2str(median(gridD.diff(gridD.diff>-inf))), ' (median of grid cells) / ',num2str(mean(gridD.diffavg(gridD.diff>-inf))),' (average)'));
    else 
        title(strcat('D = ',num2str(median(gridD.diff(gridD.diff>-inf))), ' (median of grid cells)'));
    end
    for i = 1:(length(gridD.gridx)-1)
		gx = [gridD.gridx(i) gridD.gridx(i) gridD.gridx(i+1) gridD.gridx(i+1)];
		for j = 1:(length(gridD.gridy)-1)
			gy = [gridD.gridy(j) gridD.gridy(j+1) gridD.gridy(j+1) gridD.gridy(j)];
			if (gridD.diff(i,j) > -inf)
				h = fill( gy,gx,gridD.diff(i,j), 'facealpha',.5, 'buttondownfc', {@plotGridCellDiff, gridD, i, j, par, grid.tracks});
				usedD = [usedD gridD.diff(i,j)];
			end
		end
	end

	hold off
end	
% plotGridCellDiff
% plot histogram and diffusion estimates
%
% Input
% - hObj  : object handle 
% - event : event data
% - gridD : output of getAvgGridCellDiffusion(...)/getEnsembleLSFDiffusion(...)
% - i,j	  : grid cell index
% - par	  : parameters (global variable)
% - tracks: [x y time id ...]
function plotGridCellDiff(hObj, event, gridD, i, j, par, tracks)
	subplot(2,2,1)
	title(strcat('Diffusion coefficient: ',num2str(gridD.diff(i,j))));
	t = tracks( (tracks(:,9) == i) & (tracks(:,10) == j) ,:);
	plotGridCellDiffHelper(t,par);
end	
	
% plotGridCellDiffHelper
% plot boxplot and histogram for diffusion in a grid cell (ensemble-averaged MSD) 
% (subplot(223) and subplot(224))
%
% Input
% - tracks: [x y time id loc.prec gridX gridY]
% - par	  : parameters (global variable)
function plotGridCellDiffHelper(tracks, par)
	% track step size considered for diffusion
	mint = 2;
	maxt = 4;

	if (size(tracks,1) <= maxt)
		subplot(2,2,3)
		delete(gca)
		subplot(2,2,4)
		delete(gca)
	end

	[D meanSD varSD maxn] = getMSDs(tracks,par,mint,maxt);

	if ( sum(isnan(meanSD)) < 0 || maxn > 1)
		% unweighted least squares fit
		Ols = getOLSDiffusion(tracks,par,mint,maxt,0);
		if (Ols.slo == -inf)
			subplot(2,2,3)
			delete(gca)
			subplot(2,2,4)
			delete(gca)
			return
		end
		
		% weighted least squares fit
		wOls = getOLSDiffusion(tracks,par,mint,maxt,1);
		% weighted least squares fit using pmin
		%wOlsPmin = getPminDiffusion(tracks,par);
		% robust fit
		Rls = getRobustFit(tracks,par,mint,maxt);
		% robust fit using pmin
		%RlsPmin = getRobustFit(tracks,par,mint,wOlsPmin.pmin);
		
		% plots
		subplot(2,2,3)
		boxplot(Ols.D(:,mint:maxt),'Labels',mint:maxt);
		
		hold on
		scatter(Rls.xstep-1,Rls.y)
		plot(1:maxt-mint+1,Ols.X*[Ols.off; Ols.slo],'Color','green'); % unweighted
		plot(1:maxt-mint+1,wOls.X*[wOls.off; wOls.slo],'Color','m'); % weighted
		plot(1:maxt-mint+1,Ols.X*[Rls.off; Rls.slo],'Color','red');  % robust fit
		%plot(1:wOlsPmin.pmin-mint+1,wOlsPmin.X*[wOlsPmin.off; wOlsPmin.slo],'Color','c');  % weighted + pmin
		%plot(1:wOlsPmin.pmin-mint+1,wOlsPmin.X*[RlsPmin.off; RlsPmin.slo],'Color','blue');  % robust fit + pmin
		
		ylabel('MSD');
		xlabel('dt');
		
		%title(strcat('D: OLS: ', num2str(Ols.slo/4), ' / wOLS: ', num2str(wOls.slo/4), ' / wOLS+pmin: ', num2str(wOlsPmin.slo/4), ' / RLS: ', num2str(Rls.slo/4), ' / RLS+pmin: ', num2str(RlsPmin.slo/4)));
		title(strcat('#DataPoints: ',int2str(size(tracks,1)),' / D: OLS: ', num2str(Ols.slo/4), ' / wOLS: ', num2str(wOls.slo/4), ' / RLS: ', num2str(Rls.slo/4)));
		%legend('unweighted','weighted',strcat('weighted + pmin = ',num2str(wOlsPmin.pmin)), 'robust fit', strcat('robust fit + pmin = ',num2str(wOlsPmin.pmin)), 'Location', 'NorthWest');
		legend('data','unweighted','weighted','robust fit', 'Location', 'NorthWest');
		hold off
		
		diffs = [];
		% histogram diffusion coefficients (track-based) in grid cell
		for i = unique(tracks(:,4))'
			t = tracks(tracks(:,4) == i,8)/4;
			diffs = [ diffs t(1) ];
		end
		subplot(2,2,4)
		hist(diffs)
		title('Diffusion Histogramm (selected cell, track-based)');
		
	else
		subplot(2,2,3)
		delete(gca)
		subplot(2,2,4)
		delete(gca)
	end
end	


    
% getMSDs
% get list of MSDs for a interval of step sizes
%
% Input
% - tracks	  : [x y time id ...]
% - par		  : parameters (global variable)
% - mint, maxt: time lag interval for MSD
%
% Output
% - D			 : squared displacements
% - meanSD, varSD: mean and variance of D
% - maxn		 : number of data points used for each time lag
function [D, meanSD, varSD, maxn]=getMSDs(tracks,par,mint,maxt)
	d = {};
	for (dt = mint:maxt)
		d{dt} = [];
    end

    
	for (id = unique(tracks(:,4))')
		t = tracks(tracks(:,4) == id,:);
		t(:,3) = t(:,3) - min(t(:,3)) + 1;
		
		for (dt = mint:maxt)
			sd = getSquaredDisplacements(t,dt);
            
			if ( length(sd) > 0)
				d{dt} = [d{dt} mean(sd)];
			end
		end
    end
        

	maxn = 0;
	for (dt = mint:maxt)
		maxn = max([ maxn length(d{dt}) ]);
	end

	meanSD = [];
	varSD = [];
	for (dt = mint:maxt)
		D(1:length(d{dt}),dt) = d{dt};
		meanSD = [meanSD mean(d{dt})];
		varSD = [varSD var(d{dt})];
	end
end
	
	
%getSquaredDisplacements
%
% Input
% - tracks : [x y time id ...]
% - dt	   : time lag
%
% Output
% - sd: Squared Dispacements
function sd=getSquaredDisplacements(tracks,dt)
	sd = [];
	for i = unique(tracks(:,3))'
		if ( sum(tracks(:,3)== i+dt) == 1)
			ti = tracks(tracks(:,3) == i,:);
			tj = tracks(tracks(:,3) == i+dt,:);
			sd = [sd ((tj(1) - ti(1))^2 + (tj(2) - ti(2))^2)];
		end
	end
end

	
% scale
% scales input data to interval [ min, max ]
%
% Input
% - x  : data vector
% - a,b: new min and max
%
% Output
% - res: scaled data vector
function res = scale(x,a,b)
	% unscaled avim2 has axis lim [0 772.5]
	res = (b-a) * x/772.5 + a;
end

	
% getGridIndex
% calculate grid cell index for each track
%
% Input
% - tracks: [x y ...]
% - dx, dy: size of grid cell
%
% Output
% - res
%      .tracks		  : [x y ... gridIDx gridIDy]
%      .binx, res.biny: grid coordinates
function res=getGridIndex(tracks,dx,dy)
	% grid node coordinates
	binx = min(tracks(:,1)):dx:max(tracks(:,1));
	biny = min(tracks(:,2)):dy:max(tracks(:,2));
	% grid index (bix,biy)
	[n,bix] = histc(tracks(:,1),binx);
	[n,biy] = histc(tracks(:,2),biny);
	% indices start with 1 (histc returns 0 for outside last grid cell)
	bix(bix==0) = max(bix) + 1;
	biy(biy==0) = max(biy) + 1;
	% tracks with grid cell index
	% t = [x y time track.id slope gridx gridy]
	res.tracks = [tracks bix biy];
	res.binx = binx
	res.biny = biny
end
	
% freeRoiwithin	
% find points in polygon using inpoly
% (http://www.mathworks.com/matlabcentral/fileexchange/10391-fast-points-in-polygon-test)
% much faster than inpolygon()
%
% Input
% - poly: ploygon coordinates
% - x, y: track coordinates
%
% Output
% - ind: boolean vector (1 == inside, 0 == outside)
function ind=freeRoiwithin(poly,x,y)
	%ind = inpolygon(x,y,poly(:,1),poly(:,2));
	ind = inpoly([x y],poly);
end


function pard=guidef

pard.analysismode.object=struct('String',{{'statistics','interactive track exploration','diffusion maps'}},'Style','popupmenu');
pard.analysismode.position=[1,1];
pard.analysismode.Width=3;
% ,'grid based diffusion coefficients'

pard.lentrackst.object=struct('String','Minimum length of tracks','Style','text');
pard.lentrackst.position=[2,1];
pard.lentrackst.Width=2.5;
pard.lentrackst.TooltipString=sprintf(['set this keyword to eliminate all trajectories with \n'...
            ' fewer than param.good valid positions.  This is useful \n'...
           'due to blinking noise particles in the data stream.']);
pard.lentracks.object=struct('String','2','Style','edit');
pard.lentracks.position=[2,3.5];
pard.lentracks.Width=.5;
pard.lentracks.TooltipString=pard.lentrackst.TooltipString;

pard.cutoffimmobilet.object=struct('String','Cutoff immobile log10(D(um^2/s))','Style','text');
pard.cutoffimmobilet.position=[3,1];
pard.cutoffimmobilet.Width=2.5;
pard.cutoffimmobilet.TooltipString=sprintf(['set this keyword to eliminate all trajectories with \n'...
            ' fewer than param.good valid positions.  This is useful \n'...
           'due to blinking noise particles in the data stream.']);
pard.cutoffimmobile.object=struct('String','-3','Style','edit');
pard.cutoffimmobile.position=[3,3.5];
pard.cutoffimmobile.Width=.5;
pard.cutoffimmobile.TooltipString=pard.cutoffimmobilet.TooltipString;

pard.saveD.object=struct('String','Write diffusion coefficient to localization data','Style','checkbox');
pard.saveD.position=[4,1];
pard.saveD.Width=3;


pard.plugininfo.description=sprintf('AnalyzeSPT');
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.Name='AnalyzeSPT';
end

