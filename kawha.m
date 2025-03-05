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
figure
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
% figure
% imageschs(sDEM)

sFD=FLOWobj(sDEM, 'preprocess', 'carve');
sS=STREAMobj(sFD, 'minarea',1000);
% y = smooth(sS,sS.y,'k',50,'nstribs',true); %smooth the planform by a 10 m window
% x = smooth(sS,sS.x,'k',50,'nstribs',true);
% curv = curvature(sS,x,y);
% 
% ty = smooth(strunky,strunky.y,'k',200,'nstribs',true); %smooth the planform by a 10 m window
% tx = smooth(strunky,strunky.x,'k',200,'nstribs',true);
% tcurv = curvature(strunky,tx,ty);

% figure
% imageschs(sDEM)
% hold on
% plot(sS)
% hold off
%% smooth, get curvature, do the same for trunk 

y = smooth(sS,sS.y,'k',500,'nstribs',true); %smooth the planform by a 10 m window
x = smooth(sS,sS.x,'k',500,'nstribs',true);
curv = curvature(sS,x,y);

strunky=trunk(sS);

ty = smooth(strunky,strunky.y,'k',500,'nstribs',true); %smooth the planform by a 10 m window
tx = smooth(strunky,strunky.x,'k',500,'nstribs',true);
tcurv = curvature(strunky,tx,ty);




sA=flowacc(sFD);
ixc = getnal(sS);
ixc(sS.ix) = sS.ixc;
sk = ksn(sS,sDEM,sA,.45, 10);
stk = ksn(strunky,sDEM,sA,.45, 10);

[lat,long,ixc2,elev,dist, area, smoothx, smoothy, ks, cv] = STREAMobj2XY(sS,ixc,sDEM,sS.distance, sA, x, y, sk, curv);
ixc = getnal(strunky);
ixc(strunky.ix) = strunky.ixc;
[tlat,tlong,tixc2,telev,tdist, tarea, tsmoothx, tsmoothy, tks, tcv] = STREAMobj2XY(strunky,ixc,sDEM,strunky.distance, sA, tx, ty, stk, tcurv);

ntribs_tot=length(lat(isnan(lat))) %check number of tribs
inan_tot=find(isnan(lat));

%% to get the nan index of every trib connected to the main stem
inan=[];
snapx_trib=[];
snapy_trib=[];
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
ypos_trib=NaN(length(inan), 1);
ypos_ms=NaN(length(inan), 1);
xpos_trib=NaN(length(inan), 1);
xpos_ms=NaN(length(inan), 1);

for i=1:length(inan)
    if ~isnan(i_ms(i))
        %for the trib first whats the angle
        if sum(isnan(smoothx((inan(i)-1):-1:(inan(i)-11))))==0
            poly1fit=polyfit(smoothx((inan(i)-1):-1:(inan(i)-11), 1), smoothy((inan(i)-1):-1:(inan(i)-11)), 1);
        end
        if sum(isnan(smoothx((inan(i)-1):-1:(inan(i)-11))))~=0  % do a less robust polyfit for those very short tribs
            poly1fit=polyfit(smoothx((inan(i)-1):-1:(inan(i)-3)), smoothy((inan(i)-1):-1:(inan(i)-3)), 1);
        end
        angle=atan(poly1fit(1));
        a_trib(i)=rad2deg(angle);
        %is it moving in a y pos or neg direction
        if sum(isnan(smoothx((inan(i)-1):-1:(inan(i)-11))))==0
            tribdify=smoothy(inan(i)-1)-smoothy(inan(i)-11);
            tribdifx=smoothx(inan(i)-1)-smoothx(inan(i)-11);

        end
        if sum(isnan(smoothx((inan(i)-1):-1:(inan(i)-11))))~=0  % do a less robust polyfit for those very short tribs
            tribdify=smoothy(inan(i)-1)-smoothy(inan(i)-3); 
            tribdifx=smoothx(inan(i)-1)-smoothx(inan(i)-3);        

        end
        if tribdify>0
           ypos_trib(i)=1; 
        end
        if tribdify<0
           ypos_trib(i)=0; 
        end
        if tribdifx>0
           xpos_trib(i)=1; 
        end
        if tribdifx<0
           xpos_trib(i)=0; 
        end

        %for the ms first whats the angle
        if i~=1 & i~=518
            poly1fit2=polyfit(tsmoothx((i_ms(i)+5):-1:(i_ms(i)-5)), tsmoothy((i_ms(i)+5):-1:(i_ms(i)-5)), 1);
        end
        if i==1 | i==518
            poly1fit2=polyfit(tsmoothx((i_ms(i)):-1:(i_ms(i)-5)), tsmoothy((i_ms(i)):-1:(i_ms(i)-5)), 1);
        end
        angle2=atan(poly1fit2(1));
        a_ms(i)=rad2deg(angle2);
        
        %is it moving in a y pos or neg direction
        if i~=1 & i~=518
            msdify=tsmoothy(i_ms(i)+5)-tsmoothy(i_ms(i)-5);
            msdifx=tsmoothx(i_ms(i)+5)-tsmoothx(i_ms(i)-5);
        end
        
        if i==1 | i==518
            msdify=tsmoothy(i_ms(i))-tsmoothy(i_ms(i)-5);
            msdifx=tsmoothx(i_ms(i))-tsmoothx(i_ms(i)-5);        

        end
        if msdify>0
           ypos_ms(i)=1;
        end
        if msdify<0
           ypos_ms(i)=0;
        end     
        if msdifx>0
           xpos_ms(i)=1;
        end
        if msdifx<0
           xpos_ms(i)=0;
        end  
    end
end
a_trib(ypos_trib==0 & xpos_trib>0)= a_trib(ypos_trib==0 & xpos_trib>0)+180;
a_trib(ypos_trib>0 & xpos_trib>0)= a_trib(ypos_trib>0 & xpos_trib>0)-180;
%a_trib(ypos_trib>0 & xpos_trib==0)= a_trib(ypos_trib>0 & xpos_trib==0)+360;

a_ms(ypos_ms==0 & xpos_ms>0)= a_ms(ypos_ms==0 & xpos_ms>0)+180;
a_ms(ypos_ms>0 & xpos_ms>0)= a_ms(ypos_ms>0 & xpos_ms>0)-180;
%a_ms(ypos_ms>0 & xpos_ms==0)= a_ms(ypos_ms>0 & xpos_ms==0)+360;

a_tribtoms=a_trib-a_ms;
%a_tribtoms(ypos_ms<0)=a_tribtoms(ypos_ms<0)*-1;
%a_tribtoms(ypos_trib<0)=a_tribtoms(ypos_trib<0)*-1;
a_tribtoms(a_tribtoms<180)=a_tribtoms(a_tribtoms<180)+360; %THIS IS IT
a_tribtoms(a_tribtoms>180)=a_tribtoms(a_tribtoms>180)-360;  %THIS TOO
%a_trib(ypos_trib<0)=a_trib(ypos_trib<0)*-1;
%%

figure
%imageschs(sDEM)
hold on
plot(smoothx, smoothy, 'k')
hold on
scatter(snapx_trib, snapy_trib, 60, a_tribtoms, 'filled')
%scatter(snapx_trib, snapy_trib, 20, a_trib, 'filled')
%scatter(snapx_trib, snapy_trib, 20, a_tribtoms, 'filled')

caxis([-180, 180])
colormap(ttscm('vik'))
h = colorbar;
h.Label.String = 'Angle (degrees)';
%% Get ksn at the bottom of the trib and drainage area 
ksn_tribbot=NaN(length(inan), 1);
da_trib=NaN(length(inan), 1);
curve_ms=NaN(length(inan), 1);
dist_tribbot=NaN(length(inan), 1);
%there's also a_tribtoms

for i=1:length(inan)
    if ~isnan(i_ms(i)) & i~=518
        ksn_tribbot(i)=mean(ks(inan(i)-30:inan(i)-2));
        if isnan(mean(ks(inan(i)-30 : inan(i)-2)))
            ksn_tribbot(i)=mean(ks(inan(i)-15 : inan(i)-2));
        end
        if isnan(mean(ks(inan(i)-15 : inan(i)-2)))
            ksn_tribbot(i)=mean(ks(inan(i)-5 : inan(i)-2));
        end
        da_trib(i)=area(inan(i)-1);
        curve_ms(i)=tcv(i_ms(i)-2);%mean(tcv((i_ms(i)-2):(i_ms(i)+2)));
        dist_tribbot(i)=tdist(i_ms(i)-1);
    end
end
%tomorrow must check if this actually separates tribs
%% binning by curv value
xxp=curve_ms(a_tribtoms>0 );
yyp=ksn_tribbot(a_tribtoms>0);
xxn=curve_ms(a_tribtoms<0 );
yyn=ksn_tribbot(a_tribtoms<0);
binnum=15;
curvbin=linspace(-.08, .08, binnum+1);
%curvebin=[-.15, -.04, -.025, -.015,-.005, 0, .005, .015, .025, .04, .15 ];
curvbin=[ -0.0800    -0.0480   -0.035   -0.022   -0.01   -0.004  0  0.004    0.01    0.022   0.035    0.0480       0.0800];
binnum=length(curvbin)-1;
negtbink=NaN(binnum, 1);
postbink=NaN(binnum, 1);
negtbinc=NaN(binnum, 1);
postbinc=NaN(binnum, 1);



dap=da_trib(a_tribtoms>0);
dan=da_trib(a_tribtoms<0);
angp=a_tribtoms(a_tribtoms>0);
angn=a_tribtoms(a_tribtoms<0);
negtbinda=NaN(binnum, 1);
postbinda=NaN(binnum, 1);
negtbinang=NaN(binnum, 1);
postbinang=NaN(binnum, 1);

for i=1:binnum
    
    binn=(xxn>curvbin(i) & xxn<curvbin(i+1));
    binp=(xxp>curvbin(i) & xxp<curvbin(i+1));
    if sum(binn)>4 %dont want to bin stuff with less than 5
        negtbink(i)=mean(yyn(binn), 'omitnan');
        negtbinc(i)=mean(xxn(binn), 'omitnan');        
        negtbinda(i)=mean(dan(binn), 'omitnan');        
        negtbinang(i)=mean((angn(binn)), 'omitnan');
    end
    
    if sum(binp)>4
        postbink(i)=mean(yyp(binp), 'omitnan');
        postbinc(i)=mean(xxp(binp), 'omitnan');        
        postbinda(i)=mean(dap(binp), 'omitnan');
        postbinang(i)=mean((angp(binp)), 'omitnan');
        
    end

end


%% curvature of ms vs  ksn of trib
figure
plot(curve_ms(a_tribtoms>0 & dist_tribbot<1.7e4), ksn_tribbot(a_tribtoms>0& dist_tribbot<1.7e4), 'o')

xx=curve_ms(a_tribtoms>0 );
yy=ksn_tribbot(a_tribtoms>0);
postrib=fit(postbinc(~isnan(postbink)), postbink(~isnan(postbink)), 'poly1')
figure
tiledlayout(2, 1)
nexttile
plot(curve_ms(a_tribtoms>0 ), ksn_tribbot(a_tribtoms>0), 'o', 'MarkerEdgeColor', [.5, .5, .5])
hold on
plot(postbinc, postbink, 'ks', 'MarkerFaceColor', 'r' )
plot(postrib)
xlabel('curvature')
ylabel('ksn')
title('positive tribs river right')
xlim([-.1, .1])

xx2=curve_ms(a_tribtoms<0 );
yy2=ksn_tribbot(a_tribtoms<0);
negtrib=fit(negtbinc(~isnan(negtbink)), negtbink(~isnan(negtbink)), 'poly1')
nexttile
plot(curve_ms(a_tribtoms<0 ), ksn_tribbot(a_tribtoms<0), 'o', 'MarkerEdgeColor', [.5, .5, .5])
hold on
plot(negtbinc, negtbink, 'ks', 'MarkerFaceColor', 'r' )
plot(negtrib)
xlabel('curvature')
ylabel('ksn')
title('negative tribs river left')
xlim([-.1, .1])
%% same plot but da
xx=curve_ms(a_tribtoms>0 );
yy=ksn_tribbot(a_tribtoms>0);
postrib=fit(postbinc(~isnan(postbinc)), postbinda(~isnan(postbinc)), 'poly1')
figure
tiledlayout(2, 1)
nexttile
plot(curve_ms(a_tribtoms>0 ), da_trib(a_tribtoms>0), 'o', 'MarkerEdgeColor', [.5, .5, .5])
hold on
plot(postbinc, postbinda, 'ks', 'MarkerFaceColor', 'r' )
plot(postrib)
xlabel('curvature')
ylabel('DA')
title('positive tribs river right')
xlim([-.1, .1])

%xx2=curve_ms(a_tribtoms<0 );
%yy2=ksn_tribbot(a_tribtoms<0);
negtrib=fit(negtbinc(~isnan(negtbinc)), negtbinda(~isnan(negtbinc)), 'poly1')
nexttile
plot(curve_ms(a_tribtoms<0 ), da_trib(a_tribtoms<0), 'o', 'MarkerEdgeColor', [.5, .5, .5])
hold on
plot(negtbinc, negtbinda, 'ks', 'MarkerFaceColor', 'r' )
plot(negtrib)
xlabel('curvature')
ylabel('DA')
title('negative tribs river left')
xlim([-.1, .1])
%% same plot but trib angle
absang=abs(a_tribtoms);
absang(absang>90)=absang(absang>90)-90;
xx=abs(absang(a_tribtoms>0));
yy=curve_ms(a_tribtoms>0 );
postrib=fit(xx(~isnan(yy)), yy(~isnan(yy)),  'poly1')
figure
tiledlayout(2, 1)
nexttile
plot(abs(absang(a_tribtoms>0)), curve_ms(a_tribtoms>0 ), 'o', 'MarkerEdgeColor', [.5, .5, .5])
hold on
%plot(postbinang, postbinc, 'ks', 'MarkerFaceColor', 'r' )
plot(postrib)
ylabel('curvature')
xlabel('Trib angle')
title('positive tribs river right')
ylim([-.1, .1])

xx=abs(a_tribtoms(a_tribtoms<0));
yy=curve_ms(a_tribtoms<0 );
negtrib=fit(xx(~isnan(yy)), yy(~isnan(yy)),  'poly1')
nexttile
plot(abs(absang(a_tribtoms<0)), curve_ms(a_tribtoms<0 ), 'o', 'MarkerEdgeColor', [.5, .5, .5])
hold on
%plot(negtbinang,negtbinc,  'ks', 'MarkerFaceColor', 'r' )
plot(negtrib)
ylabel('curvature')
xlabel('Trib angle')
title('negative tribs river left')
ylim([-.1, .1])

%% trib angle ksn

xx=(a_tribtoms(a_tribtoms>0));
yy=ksn_tribbot(a_tribtoms>0 );
postrib=fit(xx(~isnan(yy)), yy(~isnan(yy)),  'poly1')
figure
tiledlayout(2, 1)
nexttile
plot((a_tribtoms(a_tribtoms>0)), ksn_tribbot(a_tribtoms>0 ), 'o', 'MarkerEdgeColor', [.5, .5, .5])
hold on
%plot(postbinang, postbinc, 'ks', 'MarkerFaceColor', 'r' )
plot(postrib)
ylabel('ksn')
xlabel('Trib angle')
title('positive tribs river right')
%ylim([-.1, .1])

xx=(a_tribtoms(a_tribtoms<0));
yy=ksn_tribbot(a_tribtoms<0 );
negtrib=fit(xx(~isnan(yy)), yy(~isnan(yy)),  'poly1')
nexttile
plot((a_tribtoms(a_tribtoms<0)), ksn_tribbot(a_tribtoms<0 ), 'o', 'MarkerEdgeColor', [.5, .5, .5])
hold on
%plot(negtbinang,negtbinc,  'ks', 'MarkerFaceColor', 'r' )
plot(negtrib)
ylabel('ksn')
xlabel('Trib angle')
title('negative tribs river left')
%ylim([-.1, .1])
%%
figure
plot(dist_tribbot, ksn_tribbot, 'o')
xlabel('Dist upstream')
ylabel('ksn of trib bottom')

figure
plot(dist_tribbot, ksn_tribbot, 'o')

%%

figure
%imageschs(sDEM)
%hold on
%plot(smoothx, smoothy, 'k')
plotc(sS, sk) %curv, sk
hold on
plotc(strunky, curv) %curv, sk

%scatter(snapx_trib, snapy_trib, 60, a_tribtoms, 'filled')
%scatter(snapx_trib, snapy_trib, 50, curve_ms, 'filled')
%scatter(snapx_trib, snapy_trib, 60, a_tribtoms, 'filled')

caxis([-.1, .1])
colormap(ttscm('vik'))
h = colorbar;
h.Label.String = 'Angle (degrees)';

%%
figure
imageschs(sDEM)
%hold on
%plot(smoothx, smoothy, 'k')

plotc(sS,sk) %curv, sk
hold on
caxis([0, 80])
colormap(jet)
h = colorbar;
%h.Label.String = 'Angle (degrees)';
%plotc(strunky, tcurv) %curv, sk
% hold on
% caxis([-.03, .03])
% colormap(ttscm('vik'))
% h = colorbar;
% h.Label.String = 'Angle (degrees)';