data="F:/nz_data/airborn_lidar/kawhatau.tif"
DEM=GRIDobj(data);	
DEM.Z(DEM.Z==0)=NaN;% for some reason NaNs keep being given a 0 value
DEM = inpaintnans(DEM);

FD=FLOWobj(DEM, 'preprocess', 'carve');
S=STREAMobj(FD, 'minarea',10000);
Strunk=klargestconncomps(S, 1);
%Strunk=trunk(Strunk)
A=flowacc(FD);
%%
figure
plotdz(Strunk, DEM)

%%
c = chitransform(Strunk,A,'mn',0.55);

figure
imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false,'ticklabel','nice');
hold on
plotc(Strunk,c)
colormap(jet)
colorbar
hold off

%%
k = ksn(Strunk,DEM,A,.45, 100);
figure
histogram(k)
%%
figure
imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false,'ticklabel','nice');
hold on
plotc(Strunk,k)
colormap(jet)
colorbar
caxis([0, 100])
hold off
%% CURVATURE


y = smooth(Strunk,Strunk.y,'k',50,'nstribs',true);
x = smooth(Strunk,Strunk.x,'k',50,'nstribs',true);
curv = curvature(Strunk,x,y);
figure
histogram(curv)
%%
figure
imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false,'ticklabel','nice');
hold on
plotc(Strunk,curv)
box on
% center 0-value in colormap
%caxis([-1 1]*max(abs(caxis)))
caxis([-0.05, .05])
colormap(ttscm('vik'))
h = colorbar;
h.Label.String = 'Curvature';

hold off
%% steps
%1. clip to just the back bedrock stream
%2. get the trunk stream
%3. set trunk to NaN in a new DEM to get separate tribs 
%4. get ksn for tribs, curv for trunk and plot together
%5. for each trib get the average direction downstream
%6. separate out the right haand and left hand tribs
%7. For each trib get the curvature of the trunk where they meet (take
%bottom node, snap to trunk stream, get curvature?)
%8. For each trib get ksn of the bottom part of the stream
%9. Compare the curvature of the trunk with the ksn of the trib for left
%side and right side tribs

%% clip to just bedrock stream and make a new stream network sS
Strunky=trunk(Strunk);
i_outlet=[1862390, 5595200];
[i_outx, iout_y]=snap2stream(Strunky, i_outlet(1), i_outlet(2));
i_out_onstream=[i_outx, iout_y];
db=drainagebasins(FD, i_outx, iout_y);
sDEM=clip(DEM, db);
figure
imageschs(sDEM)

sFD=FLOWobj(sDEM, 'preprocess', 'carve');
sS=STREAMobj(sFD, 'minarea',1000);
figure
imageschs(sDEM)
hold on
plot(sS)
hold off
%% make a trunk 
strunky=trunk(sS);

sDEMnan=sDEM;

sA=flowacc(sFD);
ixc = getnal(strunky);
ixc(strunky.ix) = strunky.ixc;
[lat,long,ixc2,elev,dist, area] = STREAMobj2XY(strunky,ixc,sDEM,strunky.distance, sA);
intrunk=ismember(sDEMnan.Z, elev)
sDEMnan.Z()

