folder='D:\Data\yeast\local\160422_Ede1-mMaple_Sla2-GFP_autoON_M2\evaluation\tiffs\';

for k=1:length(se.sites)
    filename=[folder sprintf('%05d_',k) 'Ede1_Sla2GFP-asc_S' num2str(se.sites(k).ID) 'C' num2str(se.sites(k).info.cell) 'F' num2str(se.sites(k).info.filenumber) '.tif'];
    imwrite(se.sites(k).evaluation.CME2DRing.imfit.image,filename,'Compression','none');
end
