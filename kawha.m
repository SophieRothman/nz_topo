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


y = smooth(Strunk,Strunk.y,'k',200,'nstribs',true);
x = smooth(Strunk,Strunk.x,'k',200,'nstribs',true);
curv = curvature(Strunk,x,y);
figure
histogram(curv)
%%
figure
imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false,'ticklabel','nice');
hold on
%plot(x, y, '.')
plotc(Strunk,curv)
%gscatter(x, y, curv)
box on
% center 0-value in colormap
%caxis([-1 1]*max(abs(caxis)))
caxis([-0.07, .07])
colormap(ttscm('vik'))
h = colorbar;
h.Label.String = 'Curvature';

hold off
%% steps
%1. clip to just the back bedrock stream
%2. get the trunk stream
%3. set Nan index of each trib that is connected to the main stem
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
y = smooth(sS,sS.y,'k',200,'nstribs',true); %smooth the planform by a 10 m window
x = smooth(sS,sS.x,'k',200,'nstribs',true);

ty = smooth(strunky,strunky.y,'k',200,'nstribs',true); %smooth the planform by a 10 m window
tx = smooth(strunky,strunky.x,'k',200,'nstribs',true);
figure
imageschs(sDEM)
hold on
plot(sS)
hold off
%% make a trunk 
strunky=trunk(sS);


sA=flowacc(sFD);
ixc = getnal(sS);
ixc(sS.ix) = sS.ixc;
[lat,long,ixc2,elev,dist, area, smoothx, smoothy] = STREAMobj2XY(sS,ixc,sDEM,sS.distance, sA, x, y);
ixc = getnal(strunky);
ixc(strunky.ix) = strunky.ixc;
[tlat,tlong,tixc2,telev,tdist, tarea, tsmoothx, tsmoothy] = STREAMobj2XY(strunky,ixc,sDEM,strunky.distance, sA, tx, ty);

ntribs_tot=length(lat(isnan(lat))) %check number of tribs
inan_tot=find(isnan(lat));
inan=[];
snapx_trib=[];
snapy_trib=[];
%% to get the nan index of every trib connected to the main stem
count=0;
for i =2:ntribs_tot
    xout=lat(inan_tot(i)-1);
    yout=long(inan_tot(i)-1);
    [snapx, snapy]=snap2stream(strunky, xout, yout);
    dist2str=sqrt(((snapx-xout).^2)+((snapy-yout).^2));
    if dist2str<3
        count=count+1;
        inan(count, 1)=inan_tot(i);
        snapx_trib(count, 1)=snapx;
        snapy_trib(count, 1)=snapy;
        
    end
end
%% get the index in the main stream (strunky) of each junction
i_ms=NaN(length(inan), 1);
for i=1:length(inan)
    ix=find(abs(tlat-lat(inan(i)-1))<5 & abs(tlong-long(inan(i)-1))<5);
    if length(ix)>1
       distix= sqrt((tlat(ix)-lat(inan(i)-1)).^2+(tlong(ix)-long(inan(i)-1)).^2);
       mindistix=find(distix==min(distix));
       ix=ix(mindistix);
    end
    if length(ix)==1
        i_ms(i)=ix;
    end

end

%% try to get the junction angle of each trib
%first get the angle of both the trib and the mainstem for each junction
a_trib=NaN(length(inan), 1);
a_ms=NaN(length(inan), 1);

for i=1:length(inan)
    if ~isnan(i_ms(i))
        if sum(isnan(smoothx((inan(i)-11):(inan(i)-1))))==0
            poly1fit=polyfit(smoothx((inan(i)-11):(inan(i)-1)), smoothy((inan(i)-11):(inan(i)-1)), 1);
        end
        if sum(isnan(smoothx((inan(i)-11):(inan(i)-1))))~=0  % do a less robust polyfit for those very short tribs
            poly1fit=polyfit(smoothx((inan(i)-4):(inan(i)-1)), smoothy((inan(i)-4):(inan(i)-1)), 1);
        end
        angle=atan(poly1fit(1));
        a_trib(i)=rad2deg(angle);
        if i~=1 & i~=518
            poly1fit2=polyfit(tsmoothx((i_ms(i)-5):(i_ms(i)+5)), tsmoothy((i_ms(i)-5):(i_ms(i)+5)), 1);
        end
        if i==1 | i==518
            poly1fit2=polyfit(tsmoothx((i_ms(i)-5):(i_ms(i))), tsmoothy((i_ms(i)-5):(i_ms(i))), 1);
        end
            angle2=atan(poly1fit2(1));
        a_ms(i)=rad2deg(angle2);

    end
end
