data="F:/nz_data/airborn_lidar/kawhatau.tif"
DEM=GRIDobj(data);
	
DEM = inpaintnans(DEM);
FD=FLOWobj(DEM, 'preprocess', 'carve');
S=STREAMobj(FD, 'minarea',500);
A=flowacc(FD);
c = chitransform(S,A,'mn',0.45);
%%
figure
imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false,'ticklabel','nice');
hold on
plotc(S,c)
colormap(jet)
colorbar
hold off
%%
figure 
imageschs(DEM)
