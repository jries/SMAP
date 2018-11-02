function analyze_SPT(tracks,par)

% 	global ind alldat
% 	if isempty(ind)
% 		ind=1
% 	end

	% type of analysis (tracks, spl_statistics, diffusion map)
% 	mode = get(par.handles.analyzemenuespecify,'Value');
mode=par.mode;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% get tracks in region of interest
	% coordinates of the freehand roi: getPosition(pos.recinfo.hroi)
	% coordinates of the rect. roi: par.roicoords
	% 							   [x topleft, y topleft, length x, length y]
	
% 	try
% 		% get coordinates of roi
% 		coordfRoi = getPosition(pos.recinfo.hroi) / par.pixelsize;
% 		% get binary vector of tracks in roi
% 		indg = freeRoiwithin(coordfRoi,pos.tracks(:,2),pos.tracks(:,1));
% 	catch
% 		% get coordinates of roi
% 		cr = par.roicoords / par.sr_pixrec;
% 		% get binary vector of tracks in roi
% 		indg = mywithin(pos.tracks(:,2),cr([1 3]),pos.tracks(:,1),cr([2 4]));
% 	end
% 
% 	if (exist('coordfRoi') == 0) % no ROI defined -> use edge detection
% 		uiwait(warndlg('Undefined ROI! Edge Detection is used.','Undefined ROI', 'modal'));
% 		indg = edgeDetection(pos,par);
% 	else
% 		if (size(coordfRoi,1) <= 4) % most probably not a free hand ROI in this case
% 			% get coordinates of roi
% 			cr = par.roicoords / par.pixelsize;
% 			% get binary vector of tracks in roi
% 			indg = mywithin(pos.tracks(:,2),cr([1 3]),pos.tracks(:,1),cr([2 4]));
% 		end
% 	end

	% ; track will return the result 'res' (sorted by track-id)
	% ; (x) (y) (time) (track-id) (nr. of photons) (background) (loc.prec)
% 	tracks=pos.tracks(indg,1:7);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	switch mode
		
		case {1,5} % plot tracks or track-based diffusion
			% x and y = coordinates of tracks for background image
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
			
			
			if mode == 1
				% plot tracks
				figure(33)
				% plot background image
				image(y,x,avim2);
				hold on
				getD(tracks,par,2);
				hold off
				axis equal
				title(pos.info.timediff)
			else
				% plot track-based diffusion
				plotTrackDiffusion(tracks,par,x,y,avim2,-2,2);
			end
			
			
		case 2 %statistics
            cutoffimmobile=par.analv3;
			
			statt = getD(tracks,par,0);
			slo = statt.slo;
			off = statt.off;
			Ds = slo/4;
			
			figure(36)
			subplot(2,3,1)
			n = 0:0.1:5;
			h = hist(Ds,n)
			bar(n,h);
			parmhat = lognfit( Ds(Ds>0) )
			exp( parmhat(1) )
			
			alldat(ind).nD = n;
			alldat(ind).hD = h;
			
			D = median(slo)/4 %2D
			title(['D=' num2str(D)])
			xlabel('D ?m^/s')
			
			s = sqrt(median(off))
			subplot(2,3,2)
			hist((off),50)
			title(['sqrt(off)*1000=' num2str(s*1000)])
			xlabel('offset')
			
			indb = slo<=0;
			Di = slo/4;
			Di(indb) = 1e-5;
			Dlg = log10(Di);
			subplot(2,3,3)
			nlg = -4:0.1:1;
			xlabel('log10 D')
			
			indnorm = find(nlg>cutoffimmobile,1,'first');
			
			h = hist(Dlg,nlg);
			alldat(ind).nDl = nlg;alldat(ind).hDl=h;
			bar(nlg,h);
			hold off
			subplot(2,3,4);
			Dlgs = cumsum(h);
			Dlgs = Dlgs-Dlgs(indnorm);
			Dlgs = Dlgs/Dlgs(end);
			plot(nlg,Dlgs);
			xlabel('cumulative log D')
			alldat(ind).nDc = nlg;alldat(ind).hDc=Dlgs;
			subplot(2,3,5)
			nlen = 1:50;
			hlen = hist(statt.lent,nlen);
			bar(nlen,hlen);
			xlabel('tracklength')
			title(['median tracklength = ' num2str(median(statt.lent)) ' total nr of tracks ' num2str(sum(hlen))]);
			alldat(ind).nlen=nlen;alldat(ind).hlen=hlen;
			
			
		case 3 %diffusion maps
			minD = 1e-3;
			
			statt = getD(tracks,par,0);
			goodt = statt.tracks(:,5)>1e-4; % only D>0
			tr = statt.tracks(goodt,:);
			[x,indx] = sort(tr(:,1));
			y = tr(indx,2);
			D = tr(indx,5);
			
			%define grid
			dg = par.analv2; % camera pixels
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
			
			xl=get(par.handles.axes_rec,'XLim')
			yl=get(par.handles.axes_rec,'YLim')
			
			bgim=pos.palmreconstruct.palmreconstruct;
			sb=size(bgim);
			
			avim=rescaleimage(bgim,par.gamma,par.quantile);
			c=max(avim(:));
			avim(avim>c)=c;avim=avim/c;
			%     avim2=zeros([size(avim) 3]);
			%     avim2(:,:,1)=avim;avim2(:,:,2)=avim;avim2(:,:,3)=avim;
			%
			%
			minDl=par.analv3;
			maxD=par.addpar.analyze.spt(1);
			
			mask=1-2*avim;
			mask(mask<0)=0;
			figure(25)
			ind0=isnan(Dmap);
			ind0=ind0|isinf(Dmap);
			Dmap(ind0)=0;
			Dmap(Dmap<=0)=1e-5;
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
			figure(24)
			imagesc(Dldraw)
			colormap jet
			colorbar
			
			[x,y]=ind2sub(size(Dl),find(ind0));
			
			zv=0.4
			
			cm=jet(256);
			Drgb=ind2rgb(uint8(Dl*255),cm);
			
			i1=sub2ind(size(Drgb),x,y,1*ones(length(x),1));
			i2=sub2ind(size(Drgb),x,y,2*ones(length(x),1));
			i3=sub2ind(size(Drgb),x,y,3*ones(length(x),1));
			Drgb(i1)=zv;Drgb(i2)=zv;Drgb(i3)=zv;
			
			bf=(1000*par.pixelsize)/par.pixrec*dg;
			
			Db=imresize(Drgb,bf,'nearest');
			sDb=size(Db);
			
			cr=par.roicoords;
			ys=round((cr(1)-xl(1))/par.pixrec*1000)
			xs=round((cr(2)-yl(1))/par.pixrec*1000)
			Dim=zeros(sb(1),sb(2),3)+zv;
			Dim(xs+1:sDb(1)+xs,ys+1:sDb(2)+ys,:)=Db;
			
			figure(25)
			
			size(Db)
			h=imshow(Dim);
			
			%      colormap jet
			%      colorbar
			hold on
			k=   imshow(0*bgim);
			hold off
			set(k,'AlphaData',mask)
			
			
		case 4 % grid based diffusion
			
			f = figure(34)
			set(gcf, 'Position', get(0,'Screensize')); % maximize
			
			% get minD (GUI)
			minD = par.analv3;
			
			
			statt = getD(tracks,par,0);	 % calculate diffusion and plot tracks
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


	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
% edgeDetection
% find tracks within the cell using edge detection
%
% Input
% - pos: tracks and images
% - par: parameters (global variable)
%
% Output
% - ind: boolean vector (1 == track inside)
function ind=edgeDetection(pos,par)
	% x and y = coordinates of tracks for background image
	c = par.cpix;
	x = ( (c(3):c(4)+1/par.refinement+.5) + 1/par.refinement ) * par.pixelsize;
	y = ( (c(1):c(2)+1/par.refinement+.5) + 1/par.refinement ) * par.pixelsize;

	avim = pos.palmreconstruct.palmreconstruct;
	c = max( avim(:) );
	% avim( avim > c ) = c;
	avim = avim / c;
	I = zeros([size(avim) 3]);
	I(:,:,1) = avim;
	I(:,:,2) = avim;
	I(:,:,3) = avim;

	figure(90)
	% plot background image
	image(y,x,I);
	hold on

	if(size(I,3)==3)
		I=rgb2gray(I);
	end

	zoom = round(1 + 2 * par.zoom);
	h=fspecial('average',zoom);
	I=im2double(I);
	I=imfilter(I,h,'conv');

	se = strel('square',zoom - 1);
	bw = bwboundaries(imclose(I,se));

	% choose largest boundary
	s = 0;
	if iscell(bw)
		for i = 1:length(bw)
			if size(bw{i},1) > s
				b = bw{i};
				s = size(b,1);
			end
		end
	else
		b = bw
	end


	b(:,1) = scale(b(:,1),min(x),max(x));
	b(:,2) = scale(b(:,2),min(y),max(y));

	fill(b(:,2),b(:,1),'r');
	getD(pos.tracks,par,2);
	hold off

	ind = freeRoiwithin(b / par.pixelsize,pos.tracks(:,1),pos.tracks(:,2));

	
	
% getAvgGridCellDiffusion
% Compute diffusion of grid cell by averaging the diffusion coefficients of 
% tracks within the cell (weighted average)
%
% Input
%  - grid  : output of getGridIndex()
%  - par   : parameters (global variable)
%  - dx, dy: size of grid cell (default = camera pixel)
%
% Output
%  - gridD
%         .diff		   : weighted avg. diffusion
%         .gridc	   : nr of pos in cell
%         .gridx, gridy: grid coordinates
function gridD=getAvgGridCellDiffusion(grid,par,dx,dy)
	% default parameter for dx and dy are camera pixels
	if nargin < 4
		dy = par.analv2;
	end
	if nargin < 3
		dx = par.analv2;
	end
	if nargin < 2
		disp('two arguments needed: grid, par')
		return
	end

	% t = [x y time id photons background loc.prec Diff gridX gridY]
	t = grid.tracks;

	gridD.diff = repmat(-inf,length(grid.binx),length(grid.biny));
	gridD.gridc = zeros(length(grid.binx),length(grid.biny));
	for gx = unique(t(:,9))'
		for gy = unique(t(:,10))'
			% all tracks in grid cell
			tr = t( (t(:,9) == gx) & (t(:,10) == gy) ,:);
			gridD.gridc(gx,gy) = size( tr,1 );
			gridD.diff(gx,gy) = 10 ^ median( log10(tr(:,8) / 4) ) ;
            gridD.diffavg(gx,gy) = 10 ^ mean( log10(tr(:,8) / 4) ) ;
		end
	end
	gridD.diff(isnan(gridD.diff)) = -inf;
	gridD.gridx = [grid.binx grid.binx(end)+dx];
	gridD.gridy = [grid.biny grid.biny(end)+dy];

	
	
% getEnsembleLSFDiffusion
% compute grid cell diffusion by least square fit of all MSDs in the grid cell 
% (ensemble MSD)
%
% Input
% - grid  : output of getGridIndex()
% - par   : parameters (global variable)
% - dx, dy: size of grid cell (default = camera pixel)
% - robust: boolean, 1= robust LS, 0 = OLS (default = 0)
%
% Output
% - gridD
%        .tracks	   : [x y time id photons background loc.prec Diff gridX gridY]
%        .diff  	   : diffusion coefficients
%        .gridx, .gridy: grid cell coordinates
function gridD=getEnsembleLSFDiffusion(grid,par,dx,dy,robust)
	if nargin < 5
		robust = 0;
	end
	% default parameter for dx and dy are camera pixels
	if nargin < 4
		dy = par.analv2;
	end
	if nargin < 3
		dx = par.analv2;
	end
	if nargin < 2
		disp('two arguments needed: grid, par')
		return
	end

	% t = [x y time id photons background loc.prec Diff gridX gridY]
	t = grid.tracks;

	gridD.diff = repmat(-inf,length(grid.binx),length(grid.biny));
	for gx = unique(t(:,9))'
		for gy = unique(t(:,10))'
			% all tracks in grid cell
			tr = t( (t(:,9) == gx) & (t(:,10) == gy) ,:);
			if robust
				diff = getRobustFit(tr,par,2,4);
			else
				diff = getOLSDiffusion(tr, par, 2, 4, 1);
			end
			gridD.diff(gx,gy) = diff.slo / 4;
		end
	end
	
	gridD.diff(gridD.diff == NaN) = -inf;
	gridD.gridx = [grid.binx grid.binx(end)+dx];
	gridD.gridy = [grid.biny grid.biny(end)+dy];
	
	
	
	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diffusion estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% getOLSDiffusion
% estimate diffusion using ordinary least squares (weighted and unweighted)
%
% Input
% - tracks	  : [x y time id ...]
% - par		  : parameters (global)
% - mint, maxt: min. and max. considered MSD time lag
% - weighted  : boolean, 1 = weighted, 0 = unweighted, default = 1
%
% Output
% - diff
%       .slo: slope of least square fit
%       .off: offset of least square fit
%       .D: squared displacements (of getMSDs)
%       .X: matrix for least squares fit
function diff=getOLSDiffusion(tracks, par, mint, maxt, weighted)
	if nargin < 5
		weighted = 1;
	end

	diff.slo = -inf;

	if (size(tracks,1) <= maxt)
		return
	end

	[D meanSD varSD maxn] = getMSDs(tracks,par,mint,maxt);

	if ( sum(isnan(meanSD)) > 0 )
		return
	end

	t = par.timediff/1000 * (mint:maxt)';
	X = [ones(size(t)) t];
	if (weighted)
		% weighted least squares fit
		% w is inverse of variance (see Saxton 1997, Appendix)
		w = varSD.^-1;
		coeffs = lscov(X,meanSD',w);
	else
		% unweighted least squares fit
		% X  *  coeffs  =  msd
		coeffs = X \ meanSD';
	end

	diff.off = coeffs(1);
	diff.slo = coeffs(2);
	diff.D = D;
	diff.X = X;


	
% getRobustFit
% estimate diffusion using robust least squares fit
%
% Input
% - tracks	  : [x y time id ...]
% - par		  : parameters (global variable)
% - mint, maxt: min. and max. considered MSD time lag
%
% Output
% - diff
%       .slo  : slope of least square fit
%       .off  : offset of robust fit
%       .x,y  : parameters for robust fit
%       .xstep: displacements for scatter plot
function diff=getRobustFit(tracks, par, mint, maxt)
	diff.slo = -inf;

	if (size(tracks,1) <= maxt)
		return
	end

	[D meanSD varSD maxn] = getMSDs(tracks,par,mint,maxt);

	if ( sum(isnan(meanSD)) > 0 )
		return
	end

	x = [];
	xstep = [];
	y = [];
	for (dt = mint:maxt)
		t = par.timediff/1000 * dt;
		y = [y; D(:,dt)];
		x = [x; repmat(t, size( D(:,dt - mint + 1) ,1) ,1)];
		xstep = [xstep; repmat(dt, size( D(:,dt - mint + 1) ,1) ,1)];
	end

	coeffs = robustfit(x,y);

	diff.off = coeffs(1);
	diff.slo = coeffs(2);
	diff.x = x;
	diff.xstep = xstep;
	diff.y = y;
	diff.D = D;
	
	
% getD
% estimate track-based diffusion using unweighted ordinary least squares
%
% Input
% - tracks		  : [x y time id ...]
% - par			  : parameters (global variable)
% - show          : 0 => no plot 		  
%	      			1 => msd vs time 
%    				2 => plot tracks 
%
% Output
% - out			  
%      .slo		  : slo(indg);
%      .off		  : off(indg);
%      .lent	  : lent(indg);
%      .tracks	  : to;
function out=getD(tracks,par,show)
	pfound = max(tracks(:,4))
	cols = jet(pfound); 		% default colormap
	off = zeros(pfound,1);		% double vector (all zeros)
	slo = off;
	lent = off;
	indg = (slo==1);			% boolean vector (all false)
	ind1 = 1;
	st = size(tracks);
	dt = par.timediff/1000
	to = tracks;
	to(:,end+1) = 0;			% add a vector with zeros

	for k = 1:pfound
		
		% increment ind1 as long as
		%  - it is greater then the number of rows in tracks
		%  - the corrisponding row has track-id smaller k
		while( ind1 < st(1) && tracks(ind1,4) < k )
			ind1 = ind1 + 1;
		end
		
		ind2 = ind1 + 1;
		
		% increment ind2 as long as
		%  - it is greater then the number of rows in tracks
		%  - the corrisponding row has track-id equal k
		while( ind2 < st(1) && tracks(ind2,4) == k)
			ind2=ind2+1;
		end
		
		% take all rows between ind1 and ind2 (-> index always k)
		ind = ind1:ind2 - 1;
		
		% par.analv1 == min L (GUI)
		if ( length(ind) >= par.analv1 && tracks(ind(1),3) > par.framestart && tracks(ind(1),3) < par.framestop)
			indg(k) = true; 					% mark row as checked
			x = tracks(ind,1) * par.pixelsize;
			y = tracks(ind,2) * par.pixelsize;
			d = zeros( length(x)-1 ,1);			% vector with length of x
			
			% track step size considered for diffusion
			mint = 2;
			maxt = 4;
			
			% d contains the sum of mean squared distances for dt = 2 : 4 (time step size)
			for (l = mint:maxt)
				d(l) = mean( (x(1:end-l)-x(l+1:end)).^2 + (y(1:end-l)-y(l+1:end)).^2 );
			end
			
			t = dt * (mint:maxt)';
			
			% X(:,1) = 1, X(:,2) = t
			X = [ones(size(t)) t];
			
			% X  *  coeffs  =  msd
			coeffs = X \ d(mint:maxt);
			slo(k) = coeffs(2);         % slope  (k)
			off(k) = coeffs(1);         % offset (d)
			lent(k) = length(ind);		% track length
			
			to(ind,end) = coeffs(2);	% add slope to tracks
			
			if (show == 2)	%plot tracks
				plot(y,x,'Color',cols(k,:))
			end
			
			if (show == 1)						%plot msd vs time
				for (l = 1:length(x)-1)
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
	out.tracks = to;

	
	
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
	



% plotTrackDiffusion
%
% Input
% - tracks   		 : [x y time id ...]
% - par	   	 		 : parameters (global variable)
% - imgx,imgy,avim2	 : background image parameters
% - minD,maxD        : min and max diffusion coefficient
function plotTrackDiffusion(tracks,par,imgx,imgy,avim2,minD,maxD)

	pfound = max(tracks(:,4));
	ind1 = 1;
	st = size(tracks);
	dt = par.timediff/1000;

	figure(33)
	clf
	% plot background image
	image(imgy,imgx,avim2);
	
	%Min
	uicontrol('Style', 'slider',...
		'Min',-8,'Max',2,'Value',minD,...
		'Position', [80 50 500 20],...
		'Callback', {@plotTrackDiffusionHelper,tracks,par,imgx,imgy,avim2,minD,maxD,1});
	%Max	
	uicontrol('Style', 'slider',...
		'Min',-4,'Max',6,'Value',maxD,...
		'Position', [80 20 500 20],...
		'Callback', {@plotTrackDiffusionHelper,tracks,par,imgx,imgy,avim2,minD,maxD,0});
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
		while( ind1 < st(1) && tracks(ind1,4) < k )
			ind1 = ind1 + 1;
		end
		
		ind2 = ind1 + 1;
		
		% increment ind2 as long as
		%  - it is greater then the number of rows in tracks
		%  - the corrisponding row has track-id equal k
		while( ind2 < st(1) && tracks(ind2,4) == k)
			ind2=ind2+1;
		end
		
		% take all rows between ind1 and ind2 (-> index always k)
		ind = ind1:ind2 - 1;
		
		% par.analv1 == min L (GUI)
		if ( length(ind) >= par.analv1 && tracks(ind(1),3) > par.framestart && tracks(ind(1),3) < par.framestop)
			indg(k) = true; 					% mark row as checked
			x = tracks(ind,1) * par.pixelsize;
			y = tracks(ind,2) * par.pixelsize;
			d = zeros( length(x)-1 ,1);			% vector with length of x
			
			% track step size considered for diffusion
			mint = 2;
			maxt = 4;
			
			% d contains the sum of mean squared distances for dt = 2 : 4 (time step size)
			for (l = mint:maxt)
				d(l) = mean( (x(1:end-l)-x(l+1:end)).^2 + (y(1:end-l)-y(l+1:end)).^2 );
			end
			
			t = dt * (mint:maxt)';
			
			% X(:,1) = 1, X(:,2) = t
			X = [ones(size(t)) t];
			
			% X  *  coeffs  =  msd
			coeffs = X \ d(mint:maxt);
			
			if ( log10(coeffs(2) / 4) < maxD && log10(coeffs(2) / 4) > minD )
                usedD = [usedD coeffs(2)/4];
				Min = min(coeffs(2)/4,Min);
				Max = max(coeffs(2)/4,Max);
				l = color_line(y,x,repmat(coeffs(2)/4,1,length(x)));
				set(l,'ButtonDownFcn',{@plotTrackD coeffs(2)/4 text})
			end
		end
	end

	s = (Max - 0) / 6;
	h = colorbar;
	set(h,'YTick',0:s:Max)
	set(h,'YLim',[0 60])
	set(h,'YTickLabel',{0:s:Max})
	set(h,'YTickMode','auto')
	set(gcf, 'Renderer', 'opengl')

	hold off
	axis equal
	axis([ min(tracks(:,2))*0.9 max(tracks(:,2))*1.1 min(tracks(:,1))*0.9 max(tracks(:,1))*1.1 ] * par.pixelsize)
	title(['D interval = [' num2str(minD) ',' num2str(maxD) '] - Median D = ' num2str(median(usedD))])
	
	
	
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

    plotTrackDiffusion(tracks,par,imgx,imgy,avim2,minD,maxD);
    
    
    
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




	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% old functions (not used anymore, but maybe intresting)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% getPminDiffusion
% estimate diffusion using pmin (optimal number of MSDs for LSF)
% see Michalet: Mean square displacement analysis of single-particle trajectories 
%     with localization error: Brownian motion in an isotropic medium, 
%     2010 (DOI: 10.1103/PhysRevE.82.041914)
%
% Input:
% - tracks: [x y time id loc.prec.]
% - par	  : parameters (global variable)
% Output:
% - a	: offset
% - b	: slope (= 4*Diffusion Coefficient)
% - pmin: number of MSD points used
function diff=getPminDiffusion(tracks,par)
	mint = 2;
	maxt = 4;
	dt = par.timediff/1000;
	t = par.timediff/1000 * (mint:maxt)';
	usedMaxT = [4];
	diff = [];
	% until convergence
	while 1
		if ( size(tracks,1) <= maxt)
			return
		end
		
		[D meanSD varSD maxn] = getMSDs(tracks,par,mint,maxt);
		
		if ( sum(isnan(meanSD)) > 0 )
			return
		end
		
		t = par.timediff/1000 * (mint:maxt)';
		
		% weighted least squares fit
		% w is inverse of variance (Saxton 1997, Appendix)
		X = [ones(size(t)) t];
		w = varSD.^-1;
		w( w == inf | isnan(w) ) = 1;
		wcoeffs = lscov(X,meanSD',w);
		
		a = wcoeffs(1);
		b = wcoeffs(2);
		
		x = a /(b * par.timediff / 1000);  % Eq. 20
		pmin = floor(2+2.7 * sqrt(x));     % Eq. 30
		
		diff.off = a;
		diff.slo = b;
		diff.X = X;
		diff.D = D;
		diff.pmin = maxt;
		
		if (x < 0 || pmin == maxt || mint == pmin)
			break
		end
		if ( sum(usedMaxT == pmin) > 0 )
			disp('no convergence');
			break
		end
		
		maxt = pmin
		usedMaxT = [ usedMaxT pmin ];
	end
