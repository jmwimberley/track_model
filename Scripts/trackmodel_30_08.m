%% Importing of outer, inner and origin path info
Outer = readtable('Track Outer.txt');
OuterLat = Outer.latitude;
OuterLong = Outer.longitude;
OuterAlt = Outer.altitude_m_;
Inner = readtable('Track Inner.txt');
InnerLat = Inner.latitude;
InnerLong = Inner.longitude;
InnerAlt = Inner.altitude_m_;
Origin = readtable('Track Origin.txt');
OriginLat = Origin.latitude;
OriginLong = Origin.longitude;
OriginAlt = Origin.altitude_m_;
%OriginAlt = *Altitude of Origin Point*;

%% Conversion of imported Lat and Long data to XY coordinates
[x0,y0] = wgs2utm(OriginLat,OriginLong);
[xI,yI] = wgs2utm(InnerLat,InnerLong);
[xO,yO] = wgs2utm(OuterLat,OuterLong);

%%Making all points and coordinates relative to the origin
xI = xI - x0;
yI = yI - y0;
xO = xO - x0;
yO = yO - y0;
RelInnerAlt = InnerAlt - OriginAlt;
RelOuterAlt = OuterAlt - OriginAlt;

%%Plotting original relative inner and outer paths
hold on
plot3(xI,yI,RelInnerAlt)
plot3(xO,yO,RelOuterAlt)
grid on
axis equal

%% Ensuring inner and outer splines are made up of the same no. of points
xyzI = interparc(1000,xI,yI,RelInnerAlt);
xyzO = interparc(1000,xO,yO,RelOuterAlt);

%%Creating track centreline
xyzC = (xyzI+xyzO)/2;

%%Plotting the 3 interparc splines
hold on
plot3(xyzI(:,1),xyzI(:,2),xyzI(:,3))
plot3(xyzO(:,1),xyzO(:,2),xyzO(:,3))
plot3(xyzC(:,1),xyzC(:,2),xyzC(:,3))
grid on
axis equal

%% For loop to find total distance of centreline in metres (standard X-Y)
%%As well as assigning track width
format longG
distance_m(1) = 0;
dist_change(1) = sqrt(((xyzC(1,1)-xyzC(end-1,1))^2)+((xyzC(1,2)-xyzC(end-1,2))^2)+((xyzC(1,3)-xyzC(end-1,3))^2));
track_width(1) = sqrt(((xyzO(1,1)-xyzI(1,1))^2)+((xyzO(1,2)-xyzI(1,2))^2)+((xyzO(1,3)-xyzI(1,3))^2));
for j=2:length(xyzC)
    dist_change(j) = sqrt(((xyzC(j,1)-xyzC(j-1,1))^2)+((xyzC(j,2)-xyzC(j-1,2))^2)+((xyzC(j,3)-xyzC(j-1,3))^2));
     distance_m(j) = distance_m(j-1) + dist_change(j);
    track_width(j) = sqrt(((xyzO(j,1)-xyzI(j,1))^2)+((xyzO(j,2)-xyzI(j,2))^2)+((xyzO(j,3)-xyzI(j,3))^2));
end

%% RIBBON PLOT OF INNER AND OUTER SPLINES - VISUAL TRACK MODEL
sx = [xyzO(:,1)'; xyzI(:,1)'];
sy = [xyzO(:,2)'; xyzI(:,2)'];
sz = [xyzO(:,3)'; xyzO(:,3)'];

hold on
h = surf(sx,sy,sz,[distance_m;distance_m]);
grid on
axis equal
colormap(hsv(256));
h.EdgeColor = 'none';

%%Trying to add a colormap that changes with relative altitude
lowest_alt = min([xyzO(:,3);xyzI(:,3)]);
highest_alt = max([xyzO(:,3);xyzI(:,3)]);

%%Display track length in m and then km to 4 d.p
format shortG
tracklength_m = distance_m(end);
tracklength_km = (tracklength_m)/1000;

%%For loop to find track width
for i=1:length(xyzC)-1
    track_width(i) = sqrt(((xyzO(j,1)-xyzI(j,1))^2)+((xyzO(j,2)-xyzI(j-1,2))^2)+((xyzO(j,3)-xyzI(j-1,3))^2));
end

%%Heading angle calculation
heading_angle_r(1) = atan((xyzC(1,2)-xyzC(end-1,2))/(xyzC(1,1)-xyzC(end-1,1)));
for ii=2:length(xyzC)
    heading_angle_r(ii) = atan((xyzC(ii,2)-xyzC(ii-1,2))/(xyzC(ii,1)-xyzC(ii-1,1)));
end
heading_angle_d = rad2deg(heading_angle_r);

%%Setting boundaries for integral 
for ij=1:length(xyzC)
    upper_boundary_x = xyzI(ij,1);
    upper_boundary_y = xyzI(ij,2);
    upper_boundary_z = xyzI(ij,3);
    lower_boundary_x = xyzO(ij,1);
    lower_boundary_y = xyzO(ij,2);
    lower_boundary_z = xyzO(ij,3);
end

%% Two curvature calculations for each point
%%Uses formula k = d(theta)/sqrt((dx/ds)^2+(dy/ds)^2)
trans_curvature(1) = ((heading_angle_r(1)-heading_angle_r(end-1))/((dist_change(1)-dist_change(end-1))))/sqrt(((xyzC(1,1)-xyzC(end-1,1))/(dist_change(1)-dist_change(end-1)))^2+((xyzC(1,2)-xyzC(end-1,2))/(dist_change(1)-dist_change(end-1)))^2);
sag_curvature(1) = ((heading_angle_r(1)-heading_angle_r(end-1))/((dist_change(1)-dist_change(end-1))))/sqrt(((xyzC(1,1)-xyzC(end-1,1))/(dist_change(1)-dist_change(end-1)))^2+((xyzC(1,3)-xyzC(end-1,3))/(dist_change(1)-dist_change(end-1)))^2);
for jj=2:length(xyzC)
    trans_curvature(jj) = ((heading_angle_r(jj)-heading_angle_r(jj-1))/((dist_change(jj)-dist_change(jj-1))))/sqrt(((xyzC(jj,1)-xyzC(jj-1,1))/(dist_change(jj)-dist_change(jj-1)))^2+((xyzC(jj,2)-xyzC(jj-1,2))/(dist_change(jj)-dist_change(jj-1)))^2);
    sag_curvature(jj) = ((heading_angle_r(jj)-heading_angle_r(jj-1))/((dist_change(jj)-dist_change(jj-1))))/sqrt(((xyzC(jj,1)-xyzC(jj-1,1))/(dist_change(jj)-dist_change(jj-1)))^2+((xyzC(jj,3)-xyzC(jj-1,3))/(dist_change(jj)-dist_change(jj-1)))^2);
end

%%Radius of curvature
rad_trans_curv = 1./trans_curvature;
rad_sag_curv = 1./sag_curvature;

%% Slope Angle at each point
for ijk=1:length(xyzC)-1
    slope_angle_r(ijk) = atan((xyzC(ijk+1,3)-xyzC(ijk,3))/(dist_change(ijk+1)));
end
slope_angle_r(length(xyzC)) = atan((xyzC(2,3)-xyzC(end,3))/(dist_change(2)));
slope_angle_d = rad2deg(slope_angle_r);

%% Banking Angle at each point
for ijkl=1:length(xyzC)
    banking_angle_r(ijkl) = atan((xyzO(ijkl,3)-xyzI(ijkl,3))/track_width(ijkl));
end
 banking_angle_d = rad2deg(banking_angle_r);

%% 2D PLOTS
figure
plot(distance_m,slope_angle_d,'LineWidth',2,'Color','k')
xlim([min(distance_m) max(distance_m)])
ylim([min(slope_angle_d)-5 max(slope_angle_d)+5])
title('Distance Vs Slope Angle of TrackName')
xlabel('Distance (m)')
ylabel('Slope Angle (Deg.)')

figure
plot(distance_m,track_width,'LineWidth',2,'Color','k')
xlim([min(distance_m) max(distance_m)])
ylim([min(track_width)-5 max(track_width)+5])
title('Distance Vs Track Width of TrackName')
xlabel('Distance (m)')
ylabel('Track Width (m)')

figure
plot(distance_m,xyzC(:,3),'LineWidth',2,'Color','k')
xlim([min(distance_m) max(distance_m)])
ylim([min(xyzC(:,3))-5 max(xyzC(:,3))+5])
title('Distance Vs Relative Altitude of TrackName')
xlabel('Distance (m)')
ylabel('Relative Altitude to Origin (m)')

%% TORSION at each point
%First derivatives loop
dxds(1) = (xyzC(1,1)-xyzC(end-1,1))/(dist_change(1)-dist_change(end-1));
dyds(1) = (xyzC(1,2)-xyzC(end-1,2))/(dist_change(1)-dist_change(end-1));
dzds(1) = (xyzC(1,3)-xyzC(end-1,3))/(dist_change(1)-dist_change(end-1));

for ii = 2:length(xyzC);
dxds(ii) = (xyzC(ii,1)-xyzC(ii-1,1))/(dist_change(ii)-dist_change(ii-1));
dyds(ii) = (xyzC(ii,2)-xyzC(ii-1,2))/(dist_change(ii)-dist_change(ii-1));
dzds(ii) = (xyzC(ii,3)-xyzC(ii-1,3))/(dist_change(ii)-dist_change(ii-1));
end

%Second derivatives loop
d2xds2(1) = (dxds(1)-dxds(end-1))/(dist_change(1)-dist_change(end-1));
d2yds2(1) = (dyds(1)-dyds(end-1))/(dist_change(1)-dist_change(end-1));
d2zds2(1) = (dzds(1)-dzds(end-1))/(dist_change(1)-dist_change(end-1));

for ii = 2:length(xyzC);
d2xds2(ii) = (dxds(ii)-dxds(ii-1))/(dist_change(ii)-dist_change(ii-1));
d2yds2(ii) = (dyds(ii)-dyds(ii-1))/(dist_change(ii)-dist_change(ii-1));
d2zds2(ii) = (dzds(ii)-dzds(ii-1))/(dist_change(ii)-dist_change(ii-1));
end

%Third derivatives loop
d3xds3(1) = (d2xds2(1)-d2xds2(end-1))/(dist_change(1)-dist_change(end-1));
d3yds3(1) = (d2yds2(1)-d2yds2(end-1))/(dist_change(1)-dist_change(end-1));
d3zds3(1) = (d2zds2(1)-d2zds2(end-1))/(dist_change(1)-dist_change(end-1));

for ii=2:length(xyzC)
d3xds3(ii) = (d2xds2(ii)-d2xds2(ii-1))/(dist_change(ii)-dist_change(ii-1));
d3yds3(ii) = (d2yds2(ii)-d2yds2(ii-1))/(dist_change(ii)-dist_change(ii-1));
d3zds3(ii) = (d2zds2(ii)-d2zds2(ii-1))/(dist_change(ii)-dist_change(ii-1));
end

%Compile matrices and multiply by squared curvature radius
diff_matrix = cell(length(xyzC),1);
for i=1:length(xyzC)
diff_matrix{i} = [dxds(i) dyds(i) dzds(i) ; d2xds2(i) d2yds2(i) d2zds2(i) ; d3xds3(i) d3yds3(i) d3zds3(i)];
torsion(i) = (rad_trans_curv(i)^2)*det(diff_matrix{i});
end

rad_tors = 1./torsion;

%% ROTATION & SKEW-SYMMETRIC TENSOR MATRICES
%Use matrix{index} to read in specific matrices
rot_matrix = cell(length(xyzC),1);
skew_sym_tens = cell(length(xyzC),1);
for i=1:length(xyzC)
    rot_matrix{i} = [cos(heading_angle_r(i)) -sin(heading_angle_r(i)) (cos(heading_angle_r(i)).*slope_angle_r(i))+(sin(heading_angle_r(i)).*banking_angle_r(i)) ; sin(heading_angle_r(i)) cos(heading_angle_r(i)) (sin(heading_angle_r(i)).*slope_angle_r(i))-(cos(heading_angle_r(i)).*banking_angle_r(i)) ; -slope_angle_r(i) banking_angle_r(i) 1];
    skew_sym_tens{i} = [0 -trans_curvature(i) sag_curvature(i) ; trans_curvature(i) 0 -torsion(i) ; -sag_curvature(i) torsion(i) 0];
end
