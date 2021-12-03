tic
%% Importing of outer, inner and origin path info
Outer = readtable('Assen Outer.txt'); %Insert Track Name
OuterLat = Outer.latitude;
OuterLong = Outer.longitude;
OuterAlt = Outer.altitude_m_;
Inner = readtable('Assen Inner.txt'); %Insert Track Name
InnerLat = Inner.latitude;
InnerLong = Inner.longitude;
InnerAlt = Inner.altitude_m_;
Origin = readtable('Assen Origin.txt'); %Insert Track Name
OriginLat = Origin.latitude;
OriginLong = Origin.longitude;
%OriginAlt = Origin.altitude_m_;
OriginAlt = 12.0;%Altitude of Origin Point

%% Conversion of imported Lat and Long data to XY coordinates
[x0,y0] = wgs2utm(OriginLat,OriginLong);
[xI,yI] = wgs2utm(InnerLat,InnerLong);
[xO,yO] = wgs2utm(OuterLat,OuterLong);

%%Making all points and coordinates relative to the origin and combining
xI = xI - x0;
yI = yI - y0;
xO = xO - x0;
yO = yO - y0;
RelInnerAlt = InnerAlt - OriginAlt;
RelOuterAlt = OuterAlt - OriginAlt;
InnerXYZ = [xI yI RelInnerAlt];
OuterXYZ = [xO yO RelOuterAlt];

%% Generating intial centre line using loop
if length(OuterXYZ) >= length(InnerXYZ)
    for k = 1 : length(OuterXYZ)
        distances = sqrt((xO(k) - xI) .^ 2 + (yO(k) - yI) .^ 2 + (RelOuterAlt(k) - RelInnerAlt) .^ 2);
        [minDistance, indexOfMin] = min(distances);
        minimum_distances(k) = minDistance;
        minimum_index(k) = indexOfMin;
        CentreX(k) = mean([xO(k), xI(indexOfMin)]);
        CentreY(k) = mean([yO(k), yI(indexOfMin)]);
        CentreZ(k) = mean([RelOuterAlt(k), RelInnerAlt(indexOfMin)]);
%         OuterX(k) = xO(k);
%         OuterY(k) = yO(k);
%         OuterZ(k) = RelOuterAlt(k);
%         InnerX(k) = xI(indexOfMin);
%         InnerY(k) = yI(indexOfMin);
%         InnerZ(k) = RelInnerAlt(indexOfMin);
    end
else
    for k = 1 : length(InnerXYZ)
        distances = sqrt((xI(k) - xO) .^ 2 + (yI(k) - yO) .^ 2 + (RelInnerAlt(k) - RelOuterAlt) .^ 2);
        [minDistance, indexOfMin] = min(distances);
        minimum_distances(k) = minDistance;
        minimum_index(k) = indexOfMin;
        CentreX(k) = mean([xI(k), xO(indexOfMin)]);
        CentreY(k) = mean([yI(k), yO(indexOfMin)]);
        CentreZ(k) = mean([RelInnerAlt(k), RelOuterAlt(indexOfMin)]);
%         InnerX(k) = xI(k);
%         InnerY(k) = yI(k);
%         InnerZ(k) = RelInnerAlt(k);
%         OuterX(k) = xO(indexOfMin);
%         OuterY(k) = yO(indexOfMin);
%         OuterZ(k) = RelOuterAlt(indexOfMin);
    end
end
CentreXYZ = [CentreX' CentreY' CentreZ'];

%% Initial track length calculation using centre line
format longG
for j=1:length(CentreXYZ)-1
    centre_dist_change(j) = sqrt(((CentreXYZ(j+1,1)-CentreXYZ(j,1))^2)+((CentreXYZ(j+1,2)-CentreXYZ(j,2))^2)+((CentreXYZ(j+1,3)-CentreXYZ(j,3))^2));
end

centre_distance_m(1) = 0;
for j=2:length(CentreXYZ)
     centre_distance_m(j) = centre_distance_m(j-1) + centre_dist_change(j-1);
end

%Remapping points to have more accurate centre line
InnerXYZ = interparc(round(centre_distance_m(end)/5),xI,yI,RelInnerAlt);
OuterXYZ = interparc(round(centre_distance_m(end)/5),xO,yO,RelOuterAlt);
xI = InnerXYZ(:,1);
yI = InnerXYZ(:,2);
RelInnerAlt = InnerXYZ(:,3);
xO = OuterXYZ(:,1);
yO = OuterXYZ(:,2);
RelOuterAlt = OuterXYZ(:,3);
%% Making centre line and track width measurements using distance of initial centre line

if length(OuterXYZ) >= length(InnerXYZ)
    for k = 1 : length(OuterXYZ)
        distances = sqrt((xO(k) - xI) .^ 2 + (yO(k) - yI) .^ 2 + (RelOuterAlt(k) - RelInnerAlt) .^ 2);
        [minDistance, indexOfMin] = min(distances);
        minimum_distances(k) = minDistance;
        minimum_index(k) = indexOfMin;
        CentreX(k) = mean([xO(k), xI(indexOfMin)]);
        CentreY(k) = mean([yO(k), yI(indexOfMin)]);
        CentreZ(k) = mean([RelOuterAlt(k), RelInnerAlt(indexOfMin)]);
        OuterX(k) = xO(k);
        OuterY(k) = yO(k);
        OuterZ(k) = RelOuterAlt(k);
        InnerX(k) = xI(indexOfMin);
        InnerY(k) = yI(indexOfMin);
        InnerZ(k) = RelInnerAlt(indexOfMin);
        tenX(k) = (0.1*xO(k)) + (0.9*xI(indexOfMin));
        tenY(k) = (0.1*yO(k)) + (0.9*yI(indexOfMin));
        tenZ(k) = (0.1*RelOuterAlt(k)) + (0.9*RelInnerAlt(indexOfMin));
        twentyX(k) = (0.2*xO(k)) + (0.8*xI(indexOfMin));
        twentyY(k) = (0.2*yO(k)) + (0.8*yI(indexOfMin));
        twentyZ(k) = (0.2*RelOuterAlt(k)) + (0.8*RelInnerAlt(indexOfMin));
        thirtyX(k) = (0.3*xO(k)) + (0.7*xI(indexOfMin));
        thirtyY(k) = (0.3*yO(k)) + (0.7*yI(indexOfMin));
        thirtyZ(k) = (0.3*RelOuterAlt(k)) + (0.7*RelInnerAlt(indexOfMin));
        fourtyX(k) = (0.4*xO(k)) + (0.6*xI(indexOfMin));
        fourtyY(k) = (0.4*yO(k)) + (0.6*yI(indexOfMin));
        fourtyZ(k) = (0.4*RelOuterAlt(k)) + (0.6*RelInnerAlt(indexOfMin));
        sixtyX(k) = (0.6*xO(k)) + (0.4*xI(indexOfMin));
        sixtyY(k) = (0.6*yO(k)) + (0.4*yI(indexOfMin));
        sixtyZ(k) = (0.6*RelOuterAlt(k)) + (0.4*RelInnerAlt(indexOfMin));
        seventyX(k) = (0.7*xO(k)) + (0.3*xI(indexOfMin));
        seventyY(k) = (0.7*yO(k)) + (0.3*yI(indexOfMin));
        seventyZ(k) = (0.7*RelOuterAlt(k)) + (0.3*RelInnerAlt(indexOfMin));
        eightyX(k) = (0.8*xO(k)) + (0.2*xI(indexOfMin));
        eightyY(k) = (0.8*yO(k)) + (0.2*yI(indexOfMin));
        eightyZ(k) = (0.8*RelOuterAlt(k)) + (0.2*RelInnerAlt(indexOfMin));
        ninetyX(k) = (0.9*xO(k)) + (0.1*xI(indexOfMin));
        ninetyY(k) = (0.9*yO(k)) + (0.1*yI(indexOfMin));
        ninetyZ(k) = (0.9*RelOuterAlt(k)) + (0.1*RelInnerAlt(indexOfMin));
    end
else
    for k = 1 : length(InnerXYZ)
        distances = sqrt((xI(k) - xO) .^ 2 + (yI(k) - yO) .^ 2 + (RelInnerAlt(k) - RelOuterAlt) .^ 2);
        [minDistance, indexOfMin] = min(distances);
        minimum_distances(k) = minDistance;
        minimum_index(k) = indexOfMin;
        CentreX(k) = mean([xI(k), xO(indexOfMin)]);
        CentreY(k) = mean([yI(k), yO(indexOfMin)]);
        CentreZ(k) = mean([RelInnerAlt(k), RelOuterAlt(indexOfMin)]);
        InnerX(k) = xI(k);
        InnerY(k) = yI(k);
        InnerZ(k) = RelInnerAlt(k);
        OuterX(k) = xO(indexOfMin);
        OuterY(k) = yO(indexOfMin);
        OuterZ(k) = RelOuterAlt(indexOfMin);
        tenX(k) = (0.1*xO(indexOfMin)) + (0.9*xI(k));
        tenY(k) = (0.1*yO(indexOfMin)) + (0.9*yI(k));
        tenZ(k) = (0.1*RelOuterAlt(indexOfMin)) + (0.9*RelInnerAlt(k));
        twentyX(k) = (0.2*xO(indexOfMin)) + (0.8*xI(k));
        twentyY(k) = (0.2*yO(indexOfMin)) + (0.8*yI(k));
        twentyZ(k) = (0.2*RelOuterAlt(indexOfMin)) + (0.8*RelInnerAlt(k));
        thirtyX(k) = (0.3*xO(indexOfMin)) + (0.7*xI(k));
        thirtyY(k) = (0.3*yO(indexOfMin)) + (0.7*yI(k));
        thirtyZ(k) = (0.3*RelOuterAlt(indexOfMin)) + (0.7*RelInnerAlt(k));
        fourtyX(k) = (0.4*xO(indexOfMin)) + (0.6*xI(k));
        fourtyY(k) = (0.4*yO(indexOfMin)) + (0.6*yI(k));
        fourtyZ(k) = (0.4*RelOuterAlt(indexOfMin)) + (0.6*RelInnerAlt(k));
        sixtyX(k) = (0.6*xO(indexOfMin)) + (0.4*xI(k));
        sixtyY(k) = (0.6*yO(indexOfMin)) + (0.4*yI(k));
        sixtyZ(k) = (0.6*RelOuterAlt(indexOfMin)) + (0.4*RelInnerAlt(k));
        seventyX(k) = (0.7*xO(indexOfMin)) + (0.3*xI(k));
        seventyY(k) = (0.7*yO(indexOfMin)) + (0.3*yI(k));
        seventyZ(k) = (0.7*RelOuterAlt(indexOfMin)) + (0.3*RelInnerAlt(k));
        eightyX(k) = (0.8*xO(indexOfMin)) + (0.2*xI(k));
        eightyY(k) = (0.8*yO(indexOfMin)) + (0.2*yI(k));
        eightyZ(k) = (0.8*RelOuterAlt(indexOfMin)) + (0.2*RelInnerAlt(k));
        ninetyX(k) = (0.9*xO(indexOfMin)) + (0.1*xI(k));
        ninetyY(k) = (0.9*yO(indexOfMin)) + (0.1*yI(k));
        ninetyZ(k) = (0.9*RelOuterAlt(indexOfMin)) + (0.1*RelInnerAlt(k));
    end
end

%Generated lines
CentreXYZ = [CentreX' CentreY' CentreZ'];
InnerXYZ = [InnerX' InnerY' InnerZ'];
tenXYZ = [tenX' tenY' tenZ'];
twentyXYZ = [twentyX' twentyY' twentyZ'];
thirtyXYZ = [thirtyX' thirtyY' thirtyZ'];
fourtyXYZ = [fourtyX' fourtyY' fourtyZ'];
sixtyXYZ = [sixtyX' sixtyY' sixtyZ'];
seventyXYZ = [seventyX' seventyY' seventyZ'];
eightyXYZ = [eightyX' eightyY' eightyZ'];
ninetyXYZ = [ninetyX' ninetyY' ninetyZ'];
OuterXYZ = [OuterX' OuterY' OuterZ'];

%% For loop to find total distance of new centreline in metres (standard X-Y)
%%As well as assigning track width
format longG
for j=1:length(CentreXYZ)-1
    centre_dist_change(j) = sqrt(((CentreXYZ(j+1,1)-CentreXYZ(j,1))^2)+((CentreXYZ(j+1,2)-CentreXYZ(j,2))^2)+((CentreXYZ(j+1,3)-CentreXYZ(j,3))^2));
end

centre_distance_m(1) = 0;
for j=2:length(CentreXYZ)
     centre_distance_m(j) = centre_distance_m(j-1) + centre_dist_change(j-1);
end

%Track width loop
track_width(1) = sqrt(((OuterXYZ(1,1)-InnerXYZ(1,1))^2)+((OuterXYZ(1,2)-InnerXYZ(1,2))^2)+((OuterXYZ(1,3)-InnerXYZ(1,3))^2));
for j=2:length(CentreXYZ)
    track_width(j) = sqrt(((OuterXYZ(j,1)-InnerXYZ(j,1))^2)+((OuterXYZ(j,2)-InnerXYZ(j,2))^2)+((OuterXYZ(j,3)-InnerXYZ(j,3))^2));
end

%Final track length
format shortG
tracklength_m = centre_distance_m(end);
tracklength_km = (tracklength_m)/1000;

%% Cubic Spline Interpolation on 3 generated lines to gain coefficients
cubic_spline_inner = cscvn(InnerXYZ(:,[1:end 1])');
cubic_spline_outer = cscvn(OuterXYZ(:,[1:end 1])');
cubic_spline_centre = cscvn(CentreXYZ(:,[1:end 1])');
cubic_spline_ten = cscvn(tenXYZ(:,[1:end 1])');
cubic_spline_twenty = cscvn(twentyXYZ(:,[1:end 1])');
cubic_spline_thirty = cscvn(thirtyXYZ(:,[1:end 1])');
cubic_spline_fourty = cscvn(fourtyXYZ(:,[1:end 1])');
cubic_spline_sixty = cscvn(sixtyXYZ(:,[1:end 1])');
cubic_spline_seventy = cscvn(seventyXYZ(:,[1:end 1])');
cubic_spline_eighty = cscvn(eightyXYZ(:,[1:end 1])');
cubic_spline_ninety = cscvn(ninetyXYZ(:,[1:end 1])');

% %Test plot to visualise all generated lines
figure
hold on
fnplt(cubic_spline_inner,'k',0.5)
fnplt(cubic_spline_outer,'k',0.5)
fnplt(cubic_spline_centre,'g')
fnplt(cubic_spline_ten)
fnplt(cubic_spline_twenty)
fnplt(cubic_spline_thirty)
fnplt(cubic_spline_fourty)
fnplt(cubic_spline_sixty)
fnplt(cubic_spline_seventy)
fnplt(cubic_spline_eighty)
fnplt(cubic_spline_ninety)
% 
% plot3(InnerXYZ(:,1),InnerXYZ(:,2),InnerXYZ(:,3),'o')
% plot3(OuterXYZ(:,1),OuterXYZ(:,2),OuterXYZ(:,3),'o')
% plot3(CentreXYZ(:,1),CentreXYZ(:,2),CentreXYZ(:,3),'o')
% % plot3(tenXYZ(:,1),tenXYZ(:,2),tenXYZ(:,3),'o')
% % plot3(twentyXYZ(:,1),twentyXYZ(:,2),twentyXYZ(:,3),'o')
% % plot3(thirtyXYZ(:,1),thirtyXYZ(:,2),thirtyXYZ(:,3),'o')
% % plot3(fourtyXYZ(:,1),fourtyXYZ(:,2),fourtyXYZ(:,3),'o')
% % plot3(sixtyXYZ(:,1),sixtyXYZ(:,2),sixtyXYZ(:,3),'o')
% % plot3(seventyXYZ(:,1),seventyXYZ(:,2),seventyXYZ(:,3),'o')
% % plot3(eightyXYZ(:,1),eightyXYZ(:,2),eightyXYZ(:,3),'o')
% % plot3(ninetyXYZ(:,1),ninetyXYZ(:,2),ninetyXYZ(:,3),'o')
% 
%Surf plot to plot final inner to outer
sx = [OuterX; InnerX];
sy = [OuterY; InnerY];
sz = [OuterZ; InnerZ];
figure
hold on
h = surf(sx,sy,sz,[centre_distance_m;centre_distance_m]);
grid on
axis equal
colormap(hsv(256));
title('Interlagos Circuit')

%% Combining coordinates

coordinatesX = [tenX; twentyX; thirtyX; fourtyX; CentreX; sixtyX; seventyX; eightyX; ninetyX];
coordinatesY = [tenY; twentyY; thirtyY; fourtyY; CentreY; sixtyY; seventyY; eightyY; ninetyY];
coordinatesZ = [tenZ; twentyZ; thirtyZ; fourtyZ; CentreZ; sixtyZ; seventyZ; eightyZ; ninetyZ];

%% Heading angle calculations for each path
for ii=1:length(CentreXYZ)-1
    centre_heading_angle_r(ii) = atan((CentreXYZ(ii+1,2)-CentreXYZ(ii,2))/(CentreXYZ(ii+1,1)-CentreXYZ(ii,1)));
    %Inner_heading_angle_r(ii) = atan((InnerXYZ(ii+1,2)-InnerXYZ(ii,2))/(InnerXYZ(ii+1,1)-InnerXYZ(ii,1)));
    ten_heading_angle_r(ii) = atan((tenXYZ(ii+1,2)-tenXYZ(ii,2))/(tenXYZ(ii+1,1)-tenXYZ(ii,1)));
    twenty_heading_angle_r(ii) = atan((twentyXYZ(ii+1,2)-twentyXYZ(ii,2))/(twentyXYZ(ii+1,1)-twentyXYZ(ii,1)));
    thirty_heading_angle_r(ii) = atan((thirtyXYZ(ii+1,2)-thirtyXYZ(ii,2))/(thirtyXYZ(ii+1,1)-thirtyXYZ(ii,1)));
    fourty_heading_angle_r(ii) = atan((fourtyXYZ(ii+1,2)-fourtyXYZ(ii,2))/(fourtyXYZ(ii+1,1)-fourtyXYZ(ii,1)));
    sixty_heading_angle_r(ii) = atan((sixtyXYZ(ii+1,2)-sixtyXYZ(ii,2))/(sixtyXYZ(ii+1,1)-sixtyXYZ(ii,1)));
    seventy_heading_angle_r(ii) = atan((seventyXYZ(ii+1,2)-seventyXYZ(ii,2))/(seventyXYZ(ii+1,1)-seventyXYZ(ii,1)));
    eighty_heading_angle_r(ii) = atan((eightyXYZ(ii+1,2)-eightyXYZ(ii,2))/(eightyXYZ(ii+1,1)-eightyXYZ(ii,1)));
    ninety_heading_angle_r(ii) = atan((ninetyXYZ(ii+1,2)-ninetyXYZ(ii,2))/(ninetyXYZ(ii+1,1)-ninetyXYZ(ii,1)));
    %Outer_heading_angle_r(ii) = atan((OuterXYZ(ii+1,2)-OuterXYZ(ii,2))/(OuterXYZ(ii+1,1)-OuterXYZ(ii,1)));
end

centre_heading_angle_r(length(CentreXYZ)) = atan((CentreXYZ(2,2)-CentreXYZ(end,2))/(CentreXYZ(2,1)-CentreXYZ(end,1)));
%Inner_heading_angle_r(length(CentreXYZ)) = atan((InnerXYZ(2,2)-InnerXYZ(end,2))/(InnerXYZ(2,1)-InnerXYZ(end,1)));
ten_heading_angle_r(length(CentreXYZ)) = atan((tenXYZ(2,2)-tenXYZ(end,2))/(tenXYZ(2,1)-tenXYZ(end,1)));
twenty_heading_angle_r(length(CentreXYZ)) = atan((twentyXYZ(2,2)-twentyXYZ(end,2))/(twentyXYZ(2,1)-twentyXYZ(end,1)));
thirty_heading_angle_r(length(CentreXYZ)) = atan((thirtyXYZ(2,2)-thirtyXYZ(end,2))/(thirtyXYZ(2,1)-thirtyXYZ(end,1)));
fourty_heading_angle_r(length(CentreXYZ)) = atan((fourtyXYZ(2,2)-fourtyXYZ(end,2))/(fourtyXYZ(2,1)-fourtyXYZ(end,1)));
sixty_heading_angle_r(length(CentreXYZ)) = atan((sixtyXYZ(2,2)-sixtyXYZ(end,2))/(sixtyXYZ(2,1)-sixtyXYZ(end,1)));
seventy_heading_angle_r(length(CentreXYZ)) = atan((seventyXYZ(2,2)-seventyXYZ(end,2))/(seventyXYZ(2,1)-seventyXYZ(end,1)));
eighty_heading_angle_r(length(CentreXYZ)) = atan((eightyXYZ(2,2)-eightyXYZ(end,2))/(eightyXYZ(2,1)-eightyXYZ(end,1)));
ninety_heading_angle_r(length(CentreXYZ)) = atan((ninetyXYZ(2,2)-ninetyXYZ(end,2))/(ninetyXYZ(2,1)-ninetyXYZ(end,1)));
%Outer_heading_angle_r(length(CentreXYZ)) = atan((OuterXYZ(2,2)-OuterXYZ(end,2))/(OuterXYZ(2,1)-OuterXYZ(end,1)));
centre_heading_angle_d = rad2deg(heading_angle_r);

%% Slope angle calculations for each path
for ijk=1:length(CentreXYZ)-1
centre_slope_angle_r(ijk) = atan((CentreXYZ(ijk+1,3)-CentreXYZ(ijk,3))/(sqrt(((CentreXYZ(ijk+1,1)-CentreXYZ(ijk,1))^2) + ((CentreXYZ(ijk+1,2)-CentreXYZ(ijk,2))^2))));
% %Inner_slope_angle_r(ijk) = atan((InnerXYZ(ijk+1,3)-InnerXYZ(ijk,3))/(sqrt(((InnerXYZ(1,1)-InnerXYZ(end-1,1))^2) + ((InnerXYZ(1,2)-InnerXYZ(ijk,2))^2))));
% ten_slope_angle_r(ijk) = atan((tenXYZ(ijk+1,3)-tenXYZ(ijk,3))/(sqrt(((tenXYZ(ijk+1,1)-tenXYZ(ijk,1))^2) + ((tenXYZ(ijk+1,2)-tenXYZ(ijk,2))^2))));
% twenty_slope_angle_r(ijk) = atan((twentyXYZ(ijk+1,3)-twentyXYZ(end-1,3))/(sqrt(((twentyXYZ(1,1)-twentyXYZ(end-1,1))^2) + ((twentyXYZ(1,2)-twentyXYZ(end-1,2))^2))));
% thirty_slope_angle_r(ijk) = atan((thirtyXYZ(ijk+1,3)-thirtyXYZ(end-1,3))/(sqrt(((thirtyXYZ(1,1)-thirtyXYZ(end-1,1))^2) + ((thirtyXYZ(1,2)-thirtyXYZ(end-1,2))^2))));
% fourty_slope_angle_r(ijk) = atan((fourtyXYZ(ijk+1,3)-fourtyXYZ(end-1,3))/(sqrt(((fourtyXYZ(1,1)-fourtyXYZ(end-1,1))^2) + ((fourtyXYZ(1,2)-fourtyXYZ(end-1,2))^2))));
% sixty_slope_angle_r(ijk) = atan((sixtyXYZ(ijk+1,3)-sixtyXYZ(end-1,3))/(sqrt(((sixtyXYZ(1,1)-sixtyXYZ(end-1,1))^2) + ((sixtyXYZ(1,2)-sixtyXYZ(end-1,2))^2))));
% seventy_slope_angle_r(ijk) = atan((seventyXYZ(ijk+1,3)-seventyXYZ(end-1,3))/(sqrt(((seventyXYZ(1,1)-seventyXYZ(end-1,1))^2) + ((seventyXYZ(1,2)-seventyXYZ(end-1,2))^2))));
% eighty_slope_angle_r(ijk) = atan((eightyXYZ(ijk+1,3)-eightyXYZ(end-1,3))/(sqrt(((eightyXYZ(1,1)-eightyXYZ(end-1,1))^2) + ((eightyXYZ(1,2)-eightyXYZ(end-1,2))^2))));
% ninety_slope_angle_r(ijk) = atan((ninetyXYZ(1,3)-ninetyXYZ(end-1,3))/(sqrt(((ninetyXYZ(1,1)-ninetyXYZ(end-1,1))^2) + ((ninetyXYZ(1,2)-ninetyXYZ(end-1,2))^2))));
% %Outer_slope_angle_r(ijk) = atan((OuterXYZ(1,2)-OuterXYZ(end-1,2))/(sqrt(((OuterXYZ(1,1)-OuterXYZ(end-1,1))^2) + ((OuterXYZ(1,2)-OuterXYZ(end-1,2))^2))));
end

centre_slope_angle_r(length(CentreXYZ)) = atan((CentreXYZ(2,3)-CentreXYZ(end,3))/(sqrt(((CentreXYZ(2,1)-CentreXYZ(end,1))^2) + ((CentreXYZ(2,2)-CentreXYZ(end,2))^2))));
% %Inner_slope_angle_r(length(CentreXYZ)) = atan((InnerXYZ(2,3)-InnerXYZ(end,3))/(sqrt(((InnerXYZ(2,1)-InnerXYZ(end,1))^2) + ((InnerXYZ(2,2)-InnerXYZ(end,2))^2))));
% ten_slope_angle_r(length(CentreXYZ)) = atan((tenXYZ(2,3)-tenXYZ(end,3))/(sqrt(((tenXYZ(2,1)-tenXYZ(end,1))^2) + ((tenXYZ(2,2)-tenXYZ(end,2))^2))));
% twenty_slope_angle_r(length(CentreXYZ)) = atan((twentyXYZ(2,3)-twentyXYZ(end,3))/(sqrt(((twentyXYZ(2,1)-twentyXYZ(end,1))^2) + ((twentyXYZ(2,2)-twentyXYZ(end,2))^2))));
% thirty_slope_angle_r(length(CentreXYZ)) = atan((thirtyXYZ(2,3)-thirtyXYZ(end,3))/(sqrt(((thirtyXYZ(2,1)-thirtyXYZ(end,1))^2) + ((thirtyXYZ(2,2)-thirtyXYZ(end,2))^2))));
% fourty_slope_angle_r(length(CentreXYZ)) = atan((fourtyXYZ(2,3)-fourtyXYZ(end,3))/(sqrt(((fourtyXYZ(2,1)-fourtyXYZ(end,1))^2) + ((fourtyXYZ(2,2)-fourtyXYZ(end,2))^2))));
% sixty_slope_angle_r(length(CentreXYZ)) = atan((sixtyXYZ(2,3)-sixtyXYZ(end,3))/(sqrt(((sixtyXYZ(2,1)-sixtyXYZ(end,1))^2) + ((sixtyXYZ(1,2)-sixtyXYZ(end,2))^2))));
% seventy_slope_angle_r(length(CentreXYZ)) = atan((seventyXYZ(2,2)-seventyXYZ(end,2))/(sqrt(((seventyXYZ(2,1)-seventyXYZ(end,1))^2) + ((seventyXYZ(2,2)-seventyXYZ(end,2))^2))));
% eighty_slope_angle_r(length(CentreXYZ)) = atan((eightyXYZ(2,2)-eightyXYZ(end,2))/(sqrt(((eightyXYZ(2,1)-eightyXYZ(end,1))^2) + ((eightyXYZ(2,2)-eightyXYZ(end,2))^2))));
% ninety_slope_angle_r(length(CentreXYZ)) = atan((ninetyXYZ(2,2)-ninetyXYZ(end,2))/(sqrt(((ninetyXYZ(2,1)-ninetyXYZ(end,1))^2) + ((ninetyXYZ(2,2)-ninetyXYZ(end,2))^2))));
% %Outer_slope_angle_r(length(CentreXYZ)) = atan((OuterXYZ(2,2)-OuterXYZ(end,2))/(sqrt(((OuterXYZ(2,1)-OuterXYZ(end,1))^2) + ((OuterXYZ(2,2)-OuterXYZ(end,2))^2))));

centre_slope_angle_d = rad2deg(centre_slope_angle_r);

%% Derivatives loops using cubic coefficients
%Centre line derivatives
for ii = 1:length(CentreXYZ)-1
syms x
fifty_cubic_func(ii) = (cubic_spline_centre.coefs(4*(ii)-3,1)*(x^3))+(cubic_spline_centre.coefs(4*(ii)-3,2)*(x^2))+(cubic_spline_centre.coefs(4*(ii)-3,3)*(x))+(cubic_spline_centre.coefs(4*(ii)-3,4));
fifty_dxds(ii) = vpa(subs(fifty_cubic_func(ii),x,CentreX(ii)));
fifty_dyds(ii) = vpa(subs(fifty_cubic_func(ii),x,CentreY(ii)));
fifty_dzds(ii) = vpa(subs(fifty_cubic_func(ii),x,CentreZ(ii)));
% fifty_diff_cubic_func(ii) = diff(fifty_cubic_func(ii));
% fifty_d2xds2(ii) = vpa(subs(fifty_diff_cubic_func(ii),x,CentreX(ii)));
% fifty_d2yds2(ii) = vpa(subs(fifty_diff_cubic_func(ii),x,CentreY(ii)));
% fifty_d2zds2(ii) = vpa(subs(fifty_diff_cubic_func(ii),x,CentreZ(ii)));
% fifty_d_d_cubic_func(ii) = diff(fifty_diff_cubic_func(ii));
% fifty_d3xds3(ii) = vpa(subs(fifty_d_d_cubic_func(ii),x,CentreX(ii)));
% fifty_d3yds3(ii) = vpa(subs(fifty_d_d_cubic_func(ii),x,CentreY(ii)));
% fifty_d3zds3(ii) = vpa(subs(fifty_d_d_cubic_func(ii),x,CentreZ(ii)));
end

%Other derivatives
% for ii = 1:length(CentreXYZ)-1
% syms x
% inner_cubic_func(ii) = (cubic_spline_inner.coefs(4*(ii)-3,1)*(x^3))+(cubic_spline_inner.coefs(4*(ii)-3,2)*(x^2))+(cubic_spline_inner.coefs(4*(ii)-3,3)*(x))+(cubic_spline_inner.coefs(4*(ii)-3,4));
% inner_dxds(ii) = vpa(subs(inner_cubic_func(ii),x,InnerX(ii)));
% inner_dyds(ii) = vpa(subs(inner_cubic_func(ii),x,InnerY(ii)));
% inner_dzds(ii) = vpa(subs(inner_cubic_func(ii),x,InnerZ(ii)));
% % inner_diff_cubic_func(ii) = diff(inner_cubic_func(ii));
% % inner_d2xds2(ii) = vpa(subs(inner_diff_cubic_func(ii),x,InnerX(ii)));
% % inner_d2yds2(ii) = vpa(subs(inner_diff_cubic_func(ii),x,InnerY(ii)));
% % inner_d2zds2(ii) = vpa(subs(inner_diff_cubic_func(ii),x,InnerZ(ii)));
% % inner_d_d_cubic_func(ii) = diff(inner_diff_cubic_func(ii));
% % inner_d3xds3(ii) = vpa(subs(inner_d_d_cubic_func(ii),x,InnerX(ii)));
% % inner_d3yds3(ii) = vpa(subs(inner_d_d_cubic_func(ii),x,InnerY(ii)));
% % inner_d3zds3(ii) = vpa(subs(inner_d_d_cubic_func(ii),x,InnerZ(ii)));
% end

%10 percent calcs
for ii = 1:length(CentreXYZ)-1
syms x
ten_cubic_func(ii) = (cubic_spline_ten.coefs(4*(ii)-3,1)*(x^3))+(cubic_spline_ten.coefs(4*(ii)-3,2)*(x^2))+(cubic_spline_ten.coefs(4*(ii)-3,3)*(x))+(cubic_spline_ten.coefs(4*(ii)-3,4));
ten_dxds(ii) = vpa(subs(ten_cubic_func(ii),x,tenX(ii)));
ten_dyds(ii) = vpa(subs(ten_cubic_func(ii),x,tenY(ii)));
ten_dzds(ii) = vpa(subs(ten_cubic_func(ii),x,tenZ(ii)));
% ten_diff_cubic_func(ii) = diff(ten_cubic_func(ii));
% ten_d2xds2(ii) = vpa(subs(ten_diff_cubic_func(ii),x,tenX(ii)));
% ten_d2yds2(ii) = vpa(subs(ten_diff_cubic_func(ii),x,tenY(ii)));
% ten_d2zds2(ii) = vpa(subs(ten_diff_cubic_func(ii),x,tenZ(ii)));
% ten_d_d_cubic_func(ii) = diff(ten_diff_cubic_func(ii));
% ten_d3xds3(ii) = vpa(subs(ten_d_d_cubic_func(ii),x,tenX(ii)));
% ten_d3yds3(ii) = vpa(subs(ten_d_d_cubic_func(ii),x,tenY(ii)));
% ten_d3zds3(ii) = vpa(subs(ten_d_d_cubic_func(ii),x,tenZ(ii)));
end

%20 percent calcs
for ii = 1:length(CentreXYZ)-1
syms x
twenty_cubic_func(ii) = (cubic_spline_twenty.coefs(4*(ii)-3,1)*(x^3))+(cubic_spline_twenty.coefs(4*(ii)-3,2)*(x^2))+(cubic_spline_twenty.coefs(4*(ii)-3,3)*(x))+(cubic_spline_twenty.coefs(4*(ii)-3,4));
twenty_dxds(ii) = vpa(subs(twenty_cubic_func(ii),x,twentyX(ii)));
twenty_dyds(ii) = vpa(subs(twenty_cubic_func(ii),x,twentyY(ii)));
twenty_dzds(ii) = vpa(subs(twenty_cubic_func(ii),x,twentyZ(ii)));
% twenty_diff_cubic_func(ii) = diff(twenty_cubic_func(ii));
% twenty_d2xds2(ii) = vpa(subs(twenty_diff_cubic_func(ii),x,twentyX(ii)));
% twenty_d2yds2(ii) = vpa(subs(twenty_diff_cubic_func(ii),x,twentyY(ii)));
% twenty_d2zds2(ii) = vpa(subs(twenty_diff_cubic_func(ii),x,twentyZ(ii)));
% twenty_d_d_cubic_func(ii) = diff(twenty_diff_cubic_func(ii));
% twenty_d3xds3(ii) = vpa(subs(twenty_d_d_cubic_func(ii),x,twentyX(ii)));
% twenty_d3yds3(ii) = vpa(subs(twenty_d_d_cubic_func(ii),x,twentyY(ii)));
% twenty_d3zds3(ii) = vpa(subs(twenty_d_d_cubic_func(ii),x,twentyZ(ii)));
end

%30 percent calcs
for ii = 1:length(CentreXYZ)-1
syms x
thirty_cubic_func(ii) = (cubic_spline_thirty.coefs(4*(ii)-3,1)*(x^3))+(cubic_spline_thirty.coefs(4*(ii)-3,2)*(x^2))+(cubic_spline_thirty.coefs(4*(ii)-3,3)*(x))+(cubic_spline_thirty.coefs(4*(ii)-3,4));
thirty_dxds(ii) = vpa(subs(thirty_cubic_func(ii),x,thirtyX(ii)));
thirty_dyds(ii) = vpa(subs(thirty_cubic_func(ii),x,thirtyY(ii)));
thirty_dzds(ii) = vpa(subs(thirty_cubic_func(ii),x,thirtyZ(ii)));
% thirty_diff_cubic_func(ii) = diff(thirty_cubic_func(ii));
% thirty_d2xds2(ii) = vpa(subs(thirty_diff_cubic_func(ii),x,thirtyX(ii)));
% thirty_d2yds2(ii) = vpa(subs(thirty_diff_cubic_func(ii),x,thirtyY(ii)));
% thirty_d2zds2(ii) = vpa(subs(thirty_diff_cubic_func(ii),x,thirtyZ(ii)));
% thirty_d_d_cubic_func(ii) = diff(thirty_diff_cubic_func(ii));
% thirty_d3xds3(ii) = vpa(subs(thirty_d_d_cubic_func(ii),x,thirtyX(ii)));
% thirty_d3yds3(ii) = vpa(subs(thirty_d_d_cubic_func(ii),x,thirtyY(ii)));
% thirty_d3zds3(ii) = vpa(subs(thirty_d_d_cubic_func(ii),x,thirtyZ(ii)));
end

%40 percent calcs
for ii = 1:length(CentreXYZ)-1
syms x
fourty_cubic_func(ii) = (cubic_spline_fourty.coefs(4*(ii)-3,1)*(x^3))+(cubic_spline_fourty.coefs(4*(ii)-3,2)*(x^2))+(cubic_spline_fourty.coefs(4*(ii)-3,3)*(x))+(cubic_spline_fourty.coefs(4*(ii)-3,4));
fourty_dxds(ii) = vpa(subs(fourty_cubic_func(ii),x,fourtyX(ii)));
fourty_dyds(ii) = vpa(subs(fourty_cubic_func(ii),x,fourtyY(ii)));
fourty_dzds(ii) = vpa(subs(fourty_cubic_func(ii),x,fourtyZ(ii)));
% fourty_diff_cubic_func(ii) = diff(fourty_cubic_func(ii));
% fourty_d2xds2(ii) = vpa(subs(fourty_diff_cubic_func(ii),x,fourtyX(ii)));
% fourty_d2yds2(ii) = vpa(subs(fourty_diff_cubic_func(ii),x,fourtyY(ii)));
% fourty_d2zds2(ii) = vpa(subs(fourty_diff_cubic_func(ii),x,fourtyZ(ii)));
% fourty_d_d_cubic_func(ii) = diff(fourty_diff_cubic_func(ii));
% fourty_d3xds3(ii) = vpa(subs(fourty_d_d_cubic_func(ii),x,fourtyX(ii)));
% fourty_d3yds3(ii) = vpa(subs(fourty_d_d_cubic_func(ii),x,fourtyY(ii)));
% fourty_d3zds3(ii) = vpa(subs(fourty_d_d_cubic_func(ii),x,fourtyZ(ii)));
end

%60 percent calcs
for ii = 1:length(CentreXYZ)-1
syms x
sixty_cubic_func(ii) = (cubic_spline_sixty.coefs(4*(ii)-3,1)*(x^3))+(cubic_spline_sixty.coefs(4*(ii)-3,2)*(x^2))+(cubic_spline_sixty.coefs(4*(ii)-3,3)*(x))+(cubic_spline_sixty.coefs(4*(ii)-3,4));
sixty_dxds(ii) = vpa(subs(sixty_cubic_func(ii),x,sixtyX(ii)));
sixty_dyds(ii) = vpa(subs(sixty_cubic_func(ii),x,sixtyY(ii)));
sixty_dzds(ii) = vpa(subs(sixty_cubic_func(ii),x,sixtyZ(ii)));
% sixty_diff_cubic_func(ii) = diff(sixty_cubic_func(ii));
% sixty_d2xds2(ii) = vpa(subs(sixty_diff_cubic_func(ii),x,sixtyX(ii)));
% sixty_d2yds2(ii) = vpa(subs(sixty_diff_cubic_func(ii),x,sixtyY(ii)));
% sixty_d2zds2(ii) = vpa(subs(sixty_diff_cubic_func(ii),x,sixtyZ(ii)));
% sixty_d_d_cubic_func(ii) = diff(sixty_diff_cubic_func(ii));
% sixty_d3xds3(ii) = vpa(subs(sixty_d_d_cubic_func(ii),x,sixtyX(ii)));
% sixty_d3yds3(ii) = vpa(subs(sixty_d_d_cubic_func(ii),x,sixtyY(ii)));
% sixty_d3zds3(ii) = vpa(subs(sixty_d_d_cubic_func(ii),x,sixtyZ(ii)));
end

%70 percent calcs
for ii = 1:length(CentreXYZ)-1
syms x
seventy_cubic_func(ii) = (cubic_spline_seventy.coefs(4*(ii)-3,1)*(x^3))+(cubic_spline_seventy.coefs(4*(ii)-3,2)*(x^2))+(cubic_spline_seventy.coefs(4*(ii)-3,3)*(x))+(cubic_spline_seventy.coefs(4*(ii)-3,4));
seventy_dxds(ii) = vpa(subs(seventy_cubic_func(ii),x,seventyX(ii)));
seventy_dyds(ii) = vpa(subs(seventy_cubic_func(ii),x,seventyY(ii)));
seventy_dzds(ii) = vpa(subs(seventy_cubic_func(ii),x,seventyZ(ii)));
% seventy_diff_cubic_func(ii) = diff(seventy_cubic_func(ii));
% seventy_d2xds2(ii) = vpa(subs(seventy_diff_cubic_func(ii),x,seventyX(ii)));
% seventy_d2yds2(ii) = vpa(subs(seventy_diff_cubic_func(ii),x,seventyY(ii)));
% seventy_d2zds2(ii) = vpa(subs(seventy_diff_cubic_func(ii),x,seventyZ(ii)));
% seventy_d_d_cubic_func(ii) = diff(seventy_diff_cubic_func(ii));
% seventy_d3xds3(ii) = vpa(subs(seventy_d_d_cubic_func(ii),x,seventyX(ii)));
% seventy_d3yds3(ii) = vpa(subs(seventy_d_d_cubic_func(ii),x,seventyY(ii)));
% seventy_d3zds3(ii) = vpa(subs(seventy_d_d_cubic_func(ii),x,seventyZ(ii)));
end

%80 percent calcs
for ii = 1:length(CentreXYZ)-1
syms x
eighty_cubic_func(ii) = (cubic_spline_eighty.coefs(4*(ii)-3,1)*(x^3))+(cubic_spline_eighty.coefs(4*(ii)-3,2)*(x^2))+(cubic_spline_eighty.coefs(4*(ii)-3,3)*(x))+(cubic_spline_eighty.coefs(4*(ii)-3,4));
eighty_dxds(ii) = vpa(subs(eighty_cubic_func(ii),x,eightyX(ii)));
eighty_dyds(ii) = vpa(subs(eighty_cubic_func(ii),x,eightyY(ii)));
eighty_dzds(ii) = vpa(subs(eighty_cubic_func(ii),x,eightyZ(ii)));
% eighty_diff_cubic_func(ii) = diff(eighty_cubic_func(ii));
% eighty_d2xds2(ii) = vpa(subs(eighty_diff_cubic_func(ii),x,eightyX(ii)));
% eighty_d2yds2(ii) = vpa(subs(eighty_diff_cubic_func(ii),x,eightyY(ii)));
% eighty_d2zds2(ii) = vpa(subs(eighty_diff_cubic_func(ii),x,eightyZ(ii)));
% eighty_d_d_cubic_func(ii) = diff(eighty_diff_cubic_func(ii));
% eighty_d3xds3(ii) = vpa(subs(eighty_d_d_cubic_func(ii),x,eightyX(ii)));
% eighty_d3yds3(ii) = vpa(subs(eighty_d_d_cubic_func(ii),x,eightyY(ii)));
% eighty_d3zds3(ii) = vpa(subs(eighty_d_d_cubic_func(ii),x,eightyZ(ii)));
end

%90 percent calcs
for ii = 1:length(CentreXYZ)-1
syms x
ninety_cubic_func(ii) = (cubic_spline_ninety.coefs(4*(ii)-3,1)*(x^3))+(cubic_spline_ninety.coefs(4*(ii)-3,2)*(x^2))+(cubic_spline_ninety.coefs(4*(ii)-3,3)*(x))+(cubic_spline_ninety.coefs(4*(ii)-3,4));
ninety_dxds(ii) = vpa(subs(ninety_cubic_func(ii),x,ninetyX(ii)));
ninety_dyds(ii) = vpa(subs(ninety_cubic_func(ii),x,ninetyY(ii)));
ninety_dzds(ii) = vpa(subs(ninety_cubic_func(ii),x,ninetyZ(ii)));
% ninety_diff_cubic_func(ii) = diff(ninety_cubic_func(ii));
% ninety_d2xds2(ii) = vpa(subs(ninety_diff_cubic_func(ii),x,ninetyX(ii)));
% ninety_d2yds2(ii) = vpa(subs(ninety_diff_cubic_func(ii),x,ninetyY(ii)));
% ninety_d2zds2(ii) = vpa(subs(ninety_diff_cubic_func(ii),x,ninetyZ(ii)));
% ninety_d_d_cubic_func(ii) = diff(ninety_diff_cubic_func(ii));
% ninety_d3xds3(ii) = vpa(subs(ninety_d_d_cubic_func(ii),x,ninetyX(ii)));
% ninety_d3yds3(ii) = vpa(subs(ninety_d_d_cubic_func(ii),x,ninetyY(ii)));
% ninety_d3zds3(ii) = vpa(subs(ninety_d_d_cubic_func(ii),x,ninetyZ(ii)));
end

%Outer calcs
% for ii = 1:length(OuterXYZ)-1
% syms x
% outer_cubic_func(ii) = (cubic_spline_sixty.coefs(4*(ii)-3,1)*(x^3))+(cubic_spline_sixty.coefs(4*(ii)-3,2)*(x^2))+(cubic_spline_sixty.coefs(4*(ii)-3,3)*(x))+(cubic_spline_sixty.coefs(4*(ii)-3,4));
% outer_dxds(ii) = vpa(subs(sixty_cubic_func(ii),x,sixtyX(ii)));
% outer_dyds(ii) = vpa(subs(sixty_cubic_func(ii),x,sixtyY(ii)));
% outer_dzds(ii) = vpa(subs(sixty_cubic_func(ii),x,sixtyZ(ii)));
% % outer_diff_cubic_func(ii) = diff(sixty_cubic_func(ii));
% % outer_d2xds2(ii) = vpa(subs(sixty_diff_cubic_func(ii),x,sixtyX(ii)));
% % outer_d2yds2(ii) = vpa(subs(sixty_diff_cubic_func(ii),x,sixtyY(ii)));
% % outer_d2zds2(ii) = vpa(subs(sixty_diff_cubic_func(ii),x,sixtyZ(ii)));
% % outer_d_d_cubic_func(ii) = diff(sixty_diff_cubic_func(ii));
% % outer_d3xds3(ii) = vpa(subs(sixty_d_d_cubic_func(ii),x,sixtyX(ii)));
% % outer_d3yds3(ii) = vpa(subs(sixty_d_d_cubic_func(ii),x,sixtyY(ii)));
% % outer_d3zds3(ii) = vpa(subs(sixty_d_d_cubic_func(ii),x,sixtyZ(ii)));
% end

%% Curvature calculations
%Transverse curvature calculations
centre_trans_curvature(1) = ((centre_heading_angle_r(1)-centre_heading_angle_r(end-1))/(sqrt((CentreX(1)-CentreX(end-1))^2+(CentreY(1)-CentreY(end-1))^2+(CentreZ(1)-CentreZ(end-1))^2)))/(sqrt((fifty_dxds(1))^2+(fifty_dyds(1))^2));
%Inner_trans_curvature(1) = ((Inner_heading_angle_r(1)-Inner_heading_angle_r(end-1))/(sqrt((CentreX(1)-CentreX(end-1))^2+(CentreY(1)-CentreY(end-1))^2+(CentreZ(1)-CentreZ(end-1))^2)))/(sqrt((Inner_dxds(1))^2+(Inner_dyds(1))^2));
ten_trans_curvature(1) = ((ten_heading_angle_r(1)-ten_heading_angle_r(end-1))/(sqrt((tenX(1)-tenX(end-1))^2+(tenY(1)-tenY(end-1))^2+(tenZ(1)-tenZ(end-1))^2)))/(sqrt((ten_dxds(1))^2+(ten_dyds(1))^2));
twenty_trans_curvature(1) = ((twenty_heading_angle_r(1)-twenty_heading_angle_r(end-1))/(sqrt((twentyX(1)-twentyX(end-1))^2+(twentyY(1)-twentyY(end-1))^2+(twentyZ(1)-twentyZ(end-1))^2)))/(sqrt((twenty_dxds(1))^2+(twenty_dyds(1))^2));
thirty_trans_curvature(1) = ((thirty_heading_angle_r(1)-thirty_heading_angle_r(end-1))/(sqrt((thirtyX(1)-thirtyX(end-1))^2+(thirtyY(1)-thirtyY(end-1))^2+(thirtyZ(1)-thirtyZ(end-1))^2)))/(sqrt((thirty_dxds(1))^2+(thirty_dyds(1))^2));
fourty_trans_curvature(1) = ((fourty_heading_angle_r(1)-fourty_heading_angle_r(end-1))/(sqrt((fourtyX(1)-fourtyX(end-1))^2+(fourtyY(1)-fourtyY(end-1))^2+(fourtyZ(1)-fourtyZ(end-1))^2)))/(sqrt((fourty_dxds(1))^2+(fourty_dyds(1))^2));
sixty_trans_curvature(1) = ((sixty_heading_angle_r(1)-sixty_heading_angle_r(end-1))/(sqrt((sixtyX(1)-sixtyX(end-1))^2+(sixtyY(1)-sixtyY(end-1))^2+(sixtyZ(1)-sixtyZ(end-1))^2)))/(sqrt((sixty_dxds(1))^2+(sixty_dyds(1))^2));
seventy_trans_curvature(1) = ((seventy_heading_angle_r(1)-seventy_heading_angle_r(end-1))/(sqrt((seventyX(1)-seventyX(end-1))^2+(seventyY(1)-seventyY(end-1))^2+(seventyZ(1)-seventyZ(end-1))^2)))/(sqrt((seventy_dxds(1))^2+(seventy_dyds(1))^2));
eighty_trans_curvature(1) = ((eighty_heading_angle_r(1)-eighty_heading_angle_r(end-1))/(sqrt((eightyX(1)-eightyX(end-1))^2+(eightyY(1)-eightyY(end-1))^2+(eightyZ(1)-eightyZ(end-1))^2)))/(sqrt((eighty_dxds(1))^2+(eighty_dyds(1))^2));
ninety_trans_curvature(1) = ((ninety_heading_angle_r(1)-ninety_heading_angle_r(end-1))/(sqrt((ninetyX(1)-ninetyX(end-1))^2+(ninetyY(1)-ninetyY(end-1))^2+(ninetyZ(1)-ninetyZ(end-1))^2)))/(sqrt((ninety_dxds(1))^2+(ninety_dyds(1))^2));
%Outer_trans_curvature(1) = ((Outer_heading_angle_r(1)-Outer_heading_angle_r(end-1))/(sqrt((OuterX(1)-OuterX(end-1))^2+(OuterY(1)-OuterY(end-1))^2+(OuterZ(1)-OuterZ(end-1))^2)))/(sqrt((Outer_dxds(1))^2+(Outer_dyds(1))^2));

centre_sag_curvature(1) = ((centre_slope_angle_r(1)-centre_slope_angle_r(end-1))/(sqrt((CentreX(1)-CentreX(end-1))^2+(CentreY(1)-CentreY(end-1))^2+(CentreZ(1)-CentreZ(end-1))^2)))/(sqrt((fifty_dxds(1))^2+(fifty_dzds(1))^2));
%Inner_sag_curvature(1) = ((Inner_slope_angle_r(1)-Inner_slope_angle_r(end-1))/(sqrt((InnerX(1)-InnerX(end-1))^2+(InnerY(1)-InnerY(end-1))^2+(InnerZ(1)-InnerZ(end-1))^2)))/(sqrt((Inner_dxds(1))^2+(Inner_dzds(1))^2));
ten_sag_curvature(1) = ((ten_slope_angle_r(1)-ten_slope_angle_r(end-1))/(sqrt((tenX(1)-tenX(end-1))^2+(tenY(1)-tenY(end-1))^2+(tenZ(1)-tenZ(end-1))^2)))/(sqrt((ten_dxds(1))^2+(ten_dzds(1))^2));
twenty_sag_curvature(1) = ((twenty_slope_angle_r(1)-twenty_slope_angle_r(end-1))/(sqrt((twentyX(1)-twentyX(end-1))^2+(twentyY(1)-twentyY(end-1))^2+(twentyZ(1)-twentyZ(end-1))^2)))/(sqrt((twenty_dxds(1))^2+(twenty_dzds(1))^2));
thirty_sag_curvature(1) = ((thirty_slope_angle_r(1)-thirty_slope_angle_r(end-1))/(sqrt((thirtyX(1)-thirtyX(end-1))^2+(thirtyY(1)-thirtyY(end-1))^2+(thirtyZ(1)-thirtyZ(end-1))^2)))/(sqrt((thirty_dxds(1))^2+(thirty_dzds(1))^2));
fourty_sag_curvature(1) = ((fourty_slope_angle_r(1)-fourty_slope_angle_r(end-1))/(sqrt((fourtyX(1)-fourtyX(end-1))^2+(fourtyY(1)-fourtyY(end-1))^2+(fourtyZ(1)-fourtyZ(end-1))^2)))/(sqrt((fourty_dxds(1))^2+(fourty_dzds(1))^2));
sixty_sag_curvature(1) = ((sixty_slope_angle_r(1)-sixty_slope_angle_r(end-1))/(sqrt((sixtyX(1)-sixtyX(end-1))^2+(sixtyY(1)-sixtyY(end-1))^2+(sixtyZ(1)-sixtyZ(end-1))^2)))/(sqrt((sixty_dxds(1))^2+(sixty_dzds(1))^2));
seventy_sag_curvature(1) = ((seventy_slope_angle_r(1)-seventy_slope_angle_r(end-1))/(sqrt((seventyX(1)-seventyX(end-1))^2+(seventyY(1)-seventyY(end-1))^2+(seventyZ(1)-seventyZ(end-1))^2)))/(sqrt((seventy_dxds(1))^2+(seventy_dzds(1))^2));
eighty_sag_curvature(1) = ((eighty_slope_angle_r(1)-eighty_slope_angle_r(end-1))/(sqrt((eightyX(1)-eightyX(end-1))^2+(eightyY(1)-eightyY(end-1))^2+(eightyZ(1)-eightyZ(end-1))^2)))/(sqrt((eighty_dxds(1))^2+(eighty_dzds(1))^2));
ninety_sag_curvature(1) = ((ninety_slope_angle_r(1)-ninety_slope_angle_r(end-1))/(sqrt((ninetyX(1)-ninetyX(end-1))^2+(ninetyY(1)-ninetyY(end-1))^2+(ninetyZ(1)-ninetyZ(end-1))^2)))/(sqrt((ninety_dxds(1))^2+(ninety_dzds(1))^2));
%Outer_sag_curvature(1) = ((Outer_slope_angle_r(1)-Outer_slope_angle_r(end-1))/(sqrt((OuterX(1)-OuterX(end-1))^2+(OuterY(1)-OuterY(end-1))^2+(OuterZ(1)-OuterZ(end-1))^2)))/(sqrt((Outer_dxds(1))^2+(Outer_dzds(1))^2));

%Curvature loops
for jj=1:length(CentreXYZ)-1
centre_trans_curvature(jj) = ((centre_heading_angle_r(jj+1)-centre_heading_angle_r(jj))/(sqrt((CentreX(jj+1)-CentreX(jj))^2+(CentreY(jj+1)-CentreY(jj))^2+(CentreZ(jj+1)-CentreZ(jj))^2)))/(sqrt((fifty_dxds(jj))^2+(fifty_dyds(jj))^2));
%Inner_trans_curvature(jj) = ((Inner_heading_angle_r(jj+1)-Inner_heading_angle_r(jj))/(sqrt((CentreX(jj+1)-CentreX(jj))^2+(CentreY(jj+1)-CentreY(jj))^2+(CentreZ(jj+1)-CentreZ(jj))^2)))/(sqrt((Inner_dxds(jj))^2+(Inner_dyds(jj))^2));
ten_trans_curvature(jj) = ((ten_heading_angle_r(jj+1)-ten_heading_angle_r(jj))/(sqrt((tenX(jj+1)-tenX(jj))^2+(tenY(jj+1)-tenY(jj))^2+(tenZ(jj+1)-tenZ(jj))^2)))/(sqrt((ten_dxds(jj))^2+(ten_dyds(jj))^2));
twenty_trans_curvature(jj) = ((twenty_heading_angle_r(jj+1)-twenty_heading_angle_r(jj))/(sqrt((twentyX(jj+1)-twentyX(jj))^2+(twentyY(jj+1)-twentyY(jj))^2+(twentyZ(jj+1)-twentyZ(jj))^2)))/(sqrt((twenty_dxds(jj))^2+(twenty_dyds(jj))^2));
thirty_trans_curvature(jj) = ((thirty_heading_angle_r(jj+1)-thirty_heading_angle_r(jj))/(sqrt((thirtyX(jj+1)-thirtyX(jj))^2+(thirtyY(jj+1)-thirtyY(jj))^2+(thirtyZ(jj+1)-thirtyZ(jj))^2)))/(sqrt((thirty_dxds(jj))^2+(thirty_dyds(jj))^2));
fourty_trans_curvature(jj) = ((fourty_heading_angle_r(jj+1)-fourty_heading_angle_r(jj))/(sqrt((fourtyX(jj+1)-fourtyX(jj))^2+(fourtyY(jj+1)-fourtyY(jj))^2+(fourtyZ(jj+1)-fourtyZ(jj))^2)))/(sqrt((fourty_dxds(jj))^2+(fourty_dyds(jj))^2));
sixty_trans_curvature(jj) = ((sixty_heading_angle_r(jj+1)-sixty_heading_angle_r(jj))/(sqrt((sixtyX(jj+1)-sixtyX(jj))^2+(sixtyY(jj+1)-sixtyY(jj))^2+(sixtyZ(jj+1)-sixtyZ(jj))^2)))/(sqrt((sixty_dxds(jj))^2+(sixty_dyds(jj))^2));
seventy_trans_curvature(jj) = ((seventy_heading_angle_r(jj+1)-seventy_heading_angle_r(jj))/(sqrt((seventyX(jj+1)-seventyX(jj))^2+(seventyY(jj+1)-seventyY(jj))^2+(seventyZ(jj+1)-seventyZ(jj))^2)))/(sqrt((seventy_dxds(jj))^2+(seventy_dyds(jj))^2));
eighty_trans_curvature(jj) = ((eighty_heading_angle_r(jj+1)-eighty_heading_angle_r(jj))/(sqrt((eightyX(jj+1)-eightyX(jj))^2+(eightyY(jj+1)-eightyY(jj))^2+(eightyZ(jj+1)-eightyZ(jj))^2)))/(sqrt((eighty_dxds(jj))^2+(eighty_dyds(jj))^2));
ninety_trans_curvature(jj) = ((ninety_heading_angle_r(jj+1)-ninety_heading_angle_r(jj))/(sqrt((ninetyX(jj+1)-ninetyX(jj))^2+(ninetyY(jj+1)-ninetyY(jj))^2+(ninetyZ(jj+1)-ninetyZ(jj))^2)))/(sqrt((ninety_dxds(jj))^2+(ninety_dyds(jj))^2));
%Outer_trans_curvature(jj) = ((Outer_heading_angle_r(jj+1)-Outer_heading_angle_r(jj))/(sqrt((OuterX(jj+1)-OuterX(jj))^2+(OuterY(jj+1)-OuterY(jj))^2+(OuterZ(jj+1)-OuterZ(jj))^2)))/(sqrt((Outer_dxds(jj))^2+(Outer_dyds(jj))^2));

centre_sag_curvature(jj) = ((centre_slope_angle_r(jj+1)-centre_slope_angle_r(jj))/(sqrt((CentreX(jj+1)-CentreX(jj))^2+(CentreY(jj+1)-CentreY(jj))^2+(CentreZ(jj+1)-CentreZ(jj))^2)))/(sqrt((fifty_dxds(jj))^2+(fifty_dzds(jj))^2));
%Inner_sag_curvature(jj) = ((Inner_slope_angle_r(jj+1)-Inner_slope_angle_r(jj))/(sqrt((InnerX(jj+1)-InnerX(jj))^2+(InnerY(jj+1)-InnerY(jj))^2+(InnerZ(jj+1)-InnerZ(jj))^2)))/(sqrt((Inner_dxds(1))^2+(Inner_dzds(1))^2));
ten_sag_curvature(jj) = ((ten_slope_angle_r(jj+1)-ten_slope_angle_r(jj))/(sqrt((tenX(jj+1)-tenX(jj))^2+(tenY(jj+1)-tenY(jj))^2+(tenZ(jj+1)-tenZ(jj))^2)))/(sqrt((ten_dxds(1))^2+(ten_dzds(1))^2));
twenty_sag_curvature(jj) = ((twenty_slope_angle_r(jj+1)-twenty_slope_angle_r(jj))/(sqrt((twentyX(jj+1)-twentyX(jj))^2+(twentyY(jj+1)-twentyY(jj))^2+(twentyZ(jj+1)-twentyZ(jj))^2)))/(sqrt((twenty_dxds(jj))^2+(twenty_dzds(jj))^2));
thirty_sag_curvature(jj) = ((thirty_slope_angle_r(jj+1)-thirty_slope_angle_r(jj))/(sqrt((thirtyX(jj+1)-thirtyX(jj))^2+(thirtyY(jj+1)-thirtyY(jj))^2+(thirtyZ(jj+1)-thirtyZ(jj))^2)))/(sqrt((thirty_dxds(jj))^2+(thirty_dzds(jj))^2));
fourty_sag_curvature(jj) = ((fourty_slope_angle_r(jj+1)-fourty_slope_angle_r(jj))/(sqrt((fourtyX(jj+1)-fourtyX(jj))^2+(fourtyY(jj+1)-fourtyY(jj))^2+(fourtyZ(jj+1)-fourtyZ(jj))^2)))/(sqrt((fourty_dxds(jj))^2+(fourty_dzds(jj))^2));
sixty_sag_curvature(jj) = ((sixty_slope_angle_r(jj+1)-sixty_slope_angle_r(jj))/(sqrt((sixtyX(jj+1)-sixtyX(jj))^2+(sixtyY(jj+1)-sixtyY(jj))^2+(sixtyZ(jj+1)-sixtyZ(jj))^2)))/(sqrt((sixty_dxds(jj))^2+(sixty_dzds(jj))^2));
seventy_sag_curvature(jj) = ((seventy_slope_angle_r(jj+1)-seventy_slope_angle_r(jj))/(sqrt((seventyX(jj+1)-seventyX(jj))^2+(seventyY(jj+1)-seventyY(jj))^2+(seventyZ(jj+1)-seventyZ(jj))^2)))/(sqrt((seventy_dxds(jj))^2+(seventy_dzds(jj))^2));
eighty_sag_curvature(jj) = ((eighty_slope_angle_r(jj+1)-eighty_slope_angle_r(jj))/(sqrt((eightyX(jj+1)-eightyX(jj))^2+(eightyY(jj+1)-eightyY(jj))^2+(eightyZ(jj+1)-eightyZ(jj))^2)))/(sqrt((eighty_dxds(jj))^2+(eighty_dzds(jj))^2));
ninety_sag_curvature(jj) = ((ninety_slope_angle_r(jj+1)-ninety_slope_angle_r(jj))/(sqrt((ninetyX(jj+1)-ninetyX(jj))^2+(ninetyY(jj+1)-ninetyY(jj))^2+(ninetyZ(jj+1)-ninetyZ(jj))^2)))/(sqrt((ninety_dxds(jj))^2+(ninety_dzds(jj))^2));
%Outer_sag_curvature(jj) = ((Outer_slope_angle_r(jj+1)-Outer_slope_angle_r(jj))/(sqrt((OuterX(jj+1)-OuterX(jj))^2+(OuterY(jj+1)-OuterY(jj))^2+(OuterZ(jj+1)-OuterZ(jj))^2)))/(sqrt((Outer_dxds(jj))^2+(Outer_dzds(jj))^2));
end


% %% Curvature loop for finding the start point
% coordinates = { tenXYZ twentyXYZ thirtyXYZ fourtyXYZ CentreXYZ sixtyXYZ seventyXYZ eightyXYZ ninetyXYZ };
% for i = 1 : 9
%     for k = 1 : 9
%         starting_trans_curvatures(i,k) = atan((coordinates{k}(1,2)-coordinates{i}(end-1,2))/(sqrt(((coordinates{k}(1,1)-coordinates{i}(end-1,1))^2)+((coordinates{k}(1,2)-coordinates{i}(end-1,2))^2)+((coordinates{k}(1,1)-coordinates{i}(end-1,1))^2))));
%         starting_sag_curvatures(i,k) = ;
%          starting_curvatures(i,k) = sqrt(starting_trans_curvatures(i,k)^2 + starting_sag_curvatures(i,k)^2);
%         [minCurvature, indexOfMin] = min(starting_curvatures);
%         minimum_curvatures(k) = minCurvature;
%         minimum_index(k) = indexOfMin;
%     end
% end

%% Combining curvatures
for i = 1:length(CentreXYZ)-1
centre_comb_curvature(i) = sqrt(centre_trans_curvature(i)^2 + centre_sag_curvature(i)^2);
ten_comb_curvature(i) = sqrt(ten_trans_curvature(i)^2 + ten_sag_curvature(i)^2);
twenty_comb_curvature(i) = sqrt(twenty_trans_curvature(i)^2 + twenty_sag_curvature(i)^2);
thirty_comb_curvature(i) = sqrt(thirty_trans_curvature(i)^2 + thirty_sag_curvature(i)^2);
fourty_comb_curvature(i) = sqrt(fourty_trans_curvature(i)^2 + fourty_sag_curvature(i)^2);
sixty_comb_curvature(i) = sqrt(sixty_trans_curvature(i)^2 + sixty_sag_curvature(i)^2);
seventy_comb_curvature(i) = sqrt(seventy_trans_curvature(i)^2 + seventy_sag_curvature(i)^2);
eighty_comb_curvature(i) = sqrt(eighty_trans_curvature(i)^2 + eighty_sag_curvature(i)^2);
ninety_comb_curvature(i) = sqrt(ninety_trans_curvature(i)^2 + ninety_sag_curvature(i)^2);
end

%Change contents in order to change which curvatures used
curvatures = [ten_comb_curvature; twenty_comb_curvature; thirty_comb_curvature; fourty_comb_curvature; centre_comb_curvature; sixty_comb_curvature; seventy_comb_curvature; eighty_comb_curvature; ninety_comb_curvature];

for i = 1:length(curvatures)
    min_curvatures(i) = min(curvatures(:,i));
    curv_idx(i) = find(curvatures(:,i) == min_curvatures(i));
    %[r(i),c(i),p(i)] = ind2sub(size(curvatures),curv_idx(i));
end

%% METHOD TO CALCULATE AND PLOT MINIMUM CURVATURE POINTS
min_curvX(1) = mean([coordinatesX(curv_idx(1),1), coordinatesX(curv_idx(end),1)]);
min_curvY(1) = mean([coordinatesY(curv_idx(1),1), coordinatesY(curv_idx(end),1)]);
min_curvZ(1) = mean([coordinatesZ(curv_idx(1),1), coordinatesZ(curv_idx(end),1)]);
for i = 2:length(CentreXYZ)-1
    min_curvX(i) = mean([coordinatesX(curv_idx(i),i), coordinatesX(curv_idx(i-1),i)]);
    min_curvY(i) = mean([coordinatesY(curv_idx(i),i), coordinatesY(curv_idx(i-1),i)]);
    min_curvZ(i) = mean([coordinatesZ(curv_idx(i),i), coordinatesZ(curv_idx(i-1),i)]);
end
min_curvX(end+1) = mean([coordinatesX(curv_idx(1),end-1), coordinatesX(curv_idx(end),end-1)]);
min_curvY(end+1) = mean([coordinatesY(curv_idx(1),end-1), coordinatesY(curv_idx(end),end-1)]);
min_curvZ(end+1) = mean([coordinatesZ(curv_idx(1),end-1), coordinatesZ(curv_idx(end),end-1)]);

Min_CurvXYZ = [min_curvX' min_curvY' min_curvZ'];
cubic_spline_min_curv = cscvn(Min_CurvXYZ(:,[1:end 1])');
% figure
% hold on
% fnplt(cubic_spline_inner,'k',0.25)
% fnplt(cubic_spline_outer,'k',0.25)
% fnplt(cubic_spline_centre,'g',0.25)
% fnplt(cubic_spline_min_curv,'r',0.5)
% title('Minimum curvature using just x-y curvature

%% METHOD TO CALCULATE MINIMUM DISTANCE PATH
%10 percent distances
format longG
for j=1:length(tenXYZ)-1
    ten_dist_change(j) = sqrt(((tenXYZ(j+1,1)-tenXYZ(j,1))^2)+((tenXYZ(j+1,2)-tenXYZ(j,2))^2)+((tenXYZ(j+1,3)-tenXYZ(j,3))^2));
end

ten_distance_m(1) = 0;
for j=2:length(CentreXYZ)
     ten_distance_m(j) = ten_distance_m(j-1) + ten_dist_change(j-1);
end

%20 percent distances
format longG
for j=1:length(twentyXYZ)-1
    twenty_dist_change(j) = sqrt(((twentyXYZ(j+1,1)-twentyXYZ(j,1))^2)+((twentyXYZ(j+1,2)-twentyXYZ(j,2))^2)+((twentyXYZ(j+1,3)-twentyXYZ(j,3))^2));
end

twenty_distance_m(1) = 0;
for j=2:length(twentyXYZ)
     twenty_distance_m(j) = twenty_distance_m(j-1) + twenty_dist_change(j-1);
end

%30 percent distances
format longG
for j=1:length(thirtyXYZ)-1
    thirty_dist_change(j) = sqrt(((thirtyXYZ(j+1,1)-thirtyXYZ(j,1))^2)+((thirtyXYZ(j+1,2)-thirtyXYZ(j,2))^2)+((thirtyXYZ(j+1,3)-thirtyXYZ(j,3))^2));
end

thirty_distance_m(1) = 0;
for j=2:length(thirtyXYZ)
     thirty_distance_m(j) = thirty_distance_m(j-1) + thirty_dist_change(j-1);
end

%40 percent distances
format longG
for j=1:length(fourtyXYZ)-1
    fourty_dist_change(j) = sqrt(((fourtyXYZ(j+1,1)-fourtyXYZ(j,1))^2)+((fourtyXYZ(j+1,2)-fourtyXYZ(j,2))^2)+((fourtyXYZ(j+1,3)-fourtyXYZ(j,3))^2));
end

fourty_distance_m(1) = 0;
for j=2:length(fourtyXYZ)
     fourty_distance_m(j) = fourty_distance_m(j-1) + fourty_dist_change(j-1);
end

%60 percent distances
format longG
for j=1:length(sixtyXYZ)-1
    sixty_dist_change(j) = sqrt(((sixtyXYZ(j+1,1)-sixtyXYZ(j,1))^2)+((sixtyXYZ(j+1,2)-sixtyXYZ(j,2))^2)+((sixtyXYZ(j+1,3)-sixtyXYZ(j,3))^2));
end

sixty_distance_m(1) = 0;
for j=2:length(sixtyXYZ)
     sixty_distance_m(j) = sixty_distance_m(j-1) + sixty_dist_change(j-1);
end

%70 percent distances
format longG
for j=1:length(seventyXYZ)-1
    seventy_dist_change(j) = sqrt(((seventyXYZ(j+1,1)-seventyXYZ(j,1))^2)+((seventyXYZ(j+1,2)-seventyXYZ(j,2))^2)+((seventyXYZ(j+1,3)-seventyXYZ(j,3))^2));
end

seventy_distance_m(1) = 0;
for j=2:length(seventyXYZ)
     seventy_distance_m(j) = seventy_distance_m(j-1) + seventy_dist_change(j-1);
end

%80 percent distances
format longG
for j=1:length(eightyXYZ)-1
    eighty_dist_change(j) = sqrt(((eightyXYZ(j+1,1)-eightyXYZ(j,1))^2)+((eightyXYZ(j+1,2)-eightyXYZ(j,2))^2)+((eightyXYZ(j+1,3)-eightyXYZ(j,3))^2));
end

eighty_distance_m(1) = 0;
for j=2:length(eightyXYZ)
     eighty_distance_m(j) = eighty_distance_m(j-1) + eighty_dist_change(j-1);
end

%90 percent distances
format longG
for j=1:length(ninetyXYZ)-1
    ninety_dist_change(j) = sqrt(((eightyXYZ(j+1,1)-ninetyXYZ(j,1))^2)+((ninetyXYZ(j+1,2)-ninetyXYZ(j,2))^2)+((ninetyXYZ(j+1,3)-ninetyXYZ(j,3))^2));
end

ninety_distance_m(1) = 0;
for j=2:length(ninetyXYZ)
    ninety_distance_m(j) = ninety_distance_m(j-1) + ninety_dist_change(j-1);
end

%Minimum distance formation
dist_changes = [ten_dist_change; twenty_dist_change; thirty_dist_change; fourty_dist_change; centre_dist_change; sixty_dist_change; seventy_dist_change; eighty_dist_change; ninety_dist_change];

for i = 1:length(dist_changes)
    min_dist_change(i) = min(dist_changes(:,i));
    dist_idx(i) = find(dist_changes(:,i) == min_dist_change(i));
    %[r(i),c(i),p(i)] = ind2sub(size(curvatures),curv_idx(i));
end


min_distX(1) = mean([coordinatesX(dist_idx(1),1), coordinatesX(dist_idx(end),1)]);
min_distY(1) = mean([coordinatesY(dist_idx(1),1), coordinatesY(dist_idx(end),1)]);
min_distZ(1) = mean([coordinatesZ(dist_idx(1),1), coordinatesZ(dist_idx(end),1)]);
for i = 2:length(CentreXYZ)-1
    min_distX(i) = mean([coordinatesX(dist_idx(i),i), coordinatesX(dist_idx(i-1),i)]);
    min_distY(i) = mean([coordinatesY(dist_idx(i),i), coordinatesY(dist_idx(i-1),i)]);
    min_distZ(i) = mean([coordinatesZ(dist_idx(i),i), coordinatesZ(dist_idx(i-1),i)]);
end
min_distX(end+1) = mean([coordinatesX(dist_idx(1),end-1), coordinatesX(dist_idx(end),end-1)]);
min_distY(end+1) = mean([coordinatesY(dist_idx(1),end-1), coordinatesY(dist_idx(end),end-1)]);
min_distZ(end+1) = mean([coordinatesZ(dist_idx(1),end-1), coordinatesZ(dist_idx(end),end-1)]);

Min_DistXYZ = [min_distX' min_distY' min_distZ'];
cubic_spline_min_dist = cscvn(Min_DistXYZ(:,[1:end 1])');
figure
hold on
fnplt(cubic_spline_inner,'k',0.25)
fnplt(cubic_spline_outer,'k',0.25)
fnplt(cubic_spline_centre,'g',0.25)
fnplt(cubic_spline_min_curv,'r',0.5)
fnplt(cubic_spline_min_dist,'b',0.5)