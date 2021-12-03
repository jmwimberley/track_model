tic
%% Importing of outer, inner and origin path info
Outer = readtable('Imola Outer.txt'); %Insert Track Name
OuterLat = Outer.latitude;
OuterLong = Outer.longitude;
OuterAlt = Outer.altitude_m_;
Inner = readtable('Imola Inner.txt'); %Insert Track Name
InnerLat = Inner.latitude;
InnerLong = Inner.longitude;
InnerAlt = Inner.altitude_m_;
Origin = readtable('Imola Origin.txt'); %Insert Track Name
OriginLat = Origin.latitude;
OriginLong = Origin.longitude;
%OriginAlt = Origin.altitude_m_;
OriginAlt = 45.8;%Altitude of Origin Point

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
InnerXYZ = interparc(round(centre_distance_m(end)),xI,yI,RelInnerAlt);
OuterXYZ = interparc(round(centre_distance_m(end)),xO,yO,RelOuterAlt);
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
    end
end

%Generated lines
CentreXYZ = [CentreX' CentreY' CentreZ'];
InnerXYZ = [InnerX' InnerY' InnerZ'];
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

sx = [OuterX; InnerX];
sy = [OuterY; InnerY];
sz = [OuterZ; InnerZ];
figure
hold on
h = surf(sx,sy,sz,[centre_distance_m;centre_distance_m]);
grid on
axis equal
colormap(hsv(256));
title('Imola Circuit')
h.EdgeColor = 'none';
colorbar

%% QUADRATIC PROGRAMMING PROBLEM
% .. inputs::

%param reftrack        array containing the reference track, i.e. a reference line and the according track widths
%                             to the right and to the left [x, y, z, w_tr_right, w_tr_left] (unit is meter, must be unclosed)

%param normvectors     normalized normal vectors for every point of the reference track [x_component, y_component]
%                             (unit is meter, must be unclosed!)


ref_vectorsX = CentreX(2:end)-CentreX(1:end-1);
ref_vectorsY = CentreY(2:end)-CentreY(1:end-1);
ref_vectorsZ = CentreZ(2:end)-CentreZ(1:end-1);
ref_vectors = [ref_vectorsX ; ref_vectorsY ; ref_vectorsZ];

%Normalised normal vectors (As 3D and of form (a,b,c), normal is (-b,a,0)
for i = 1:length(CentreXYZ)-1
    norm_vectorsX(i) = -ref_vectorsY(i)/sqrt(((CentreX(i+1)-CentreX(i))^2)+((CentreY(i+1)-CentreY(i))^2)+((CentreZ(i+1)-CentreZ(i))^2));
    norm_vectorsY(i) = ref_vectorsX(i)/sqrt(((CentreX(i+1)-CentreX(i))^2)+((CentreY(i+1)-CentreY(i))^2)+((CentreZ(i+1)-CentreZ(i))^2));
    norm_vectorsZ(i) = 0;
end
 norm_vectors = [norm_vectorsX ; norm_vectorsY ; norm_vectorsZ];

%% INPUTS
% ref_track - Requires the x, y and z coordinates of the reference (in this case centre) line, to then build a double that takes the form uses [x, y, z, w_tr_right, w_tr_left] in metres
%            The other requirement is arrays for w_tr_right and w_tr_left,
%            which uses the track width calculation of the inner and outer
%            paths and gets added to either side of the centre line in
%            order to map the varying track width.
% normvectors - Normalised normal vectors for every point of the reference
%               track [x_component, y_component, z_component] in metres.

%% Building reftrack

w_tr_right = track_width / 2;
w_tr_left = track_width / 2;
w_tr_right(end) = [];
w_tr_left(end) = [];

ref_track = [ref_vectorsX ; ref_vectorsY ; ref_vectorsZ ; w_tr_right ; w_tr_left];

%% Setting up final matrices for the solver
no_points = length(ref_track);
H = zeros(no_points);
f = zeros(1,no_points);

H(1,1) = 4*(norm_vectorsX(1)^2 + norm_vectorsY(1)^2 + norm_vectorsZ(1)^2);
H(1,2) = 0.5*2*(-2*norm_vectorsX(1)*norm_vectorsX(2)-(2*norm_vectorsY(1)*norm_vectorsY(2))-(2*norm_vectorsZ(1)*norm_vectorsZ(2)));
H(2,1) = H(1,2);
f(1) = -2*norm_vectorsX(1)*ref_track(1,end)-2*norm_vectorsY(1)*ref_track(2,end)-2*norm_vectorsZ(1)*ref_track(3,end)+4*norm_vectorsX(1)*ref_track(1,1)-2*norm_vectorsX(1)*ref_track(1,2)+4*norm_vectorsY(1)*ref_track(2,1)-2*norm_vectorsY(1)*ref_track(2,2)+4*norm_vectorsZ(1)*ref_track(3,1)-2*norm_vectorsZ(1)*ref_track(3,2);
for i = 2:no_points-1
   H(i,i) = 4*(norm_vectorsX(i)^2 + norm_vectorsY(i)^2 + norm_vectorsZ(i)^2);
   H(i,i+1) = 0.5*2*(-2*norm_vectorsX(i)*norm_vectorsX(i+1)-(2*norm_vectorsY(i)*norm_vectorsY(i+1))-(2*norm_vectorsZ(i)*norm_vectorsZ(i+1)));
   H(i+1,i) = H(i,i+1);
   f(i)= -2*norm_vectorsX(i)*ref_track(1,i-1)-2*norm_vectorsY(i)*ref_track(2,i-1)-2*norm_vectorsZ(i)*ref_track(3,i-1)+4*norm_vectorsX(i)*ref_track(1,i)-2*norm_vectorsX(i)*ref_track(1,i+1)+4*norm_vectorsY(i)*ref_track(2,i)-2*norm_vectorsY(i)*ref_track(2,i+1)+4*norm_vectorsZ(i)*ref_track(3,i)-2*norm_vectorsZ(i)*ref_track(3,i+1);  
end
H(no_points,no_points) = 4*(norm_vectorsX(no_points)^2 + norm_vectorsY(no_points)^2 + norm_vectorsZ(no_points)^2);
H(no_points,1) = 0.5*2*(-2*norm_vectorsX(no_points)*norm_vectorsX(1)-(2*norm_vectorsY(no_points)*norm_vectorsY(1))-(2*norm_vectorsZ(no_points)*norm_vectorsZ(1)));
H(1,no_points) = H(no_points,1);
f(no_points) = -2*norm_vectorsX(no_points)*ref_track(1,no_points-1)-2*norm_vectorsY(no_points)*ref_track(2,no_points-1)-2*norm_vectorsZ(no_points)*ref_track(3,no_points-1)+4*norm_vectorsX(no_points)*ref_track(1,no_points)+4*norm_vectorsY(no_points)*ref_track(2,no_points)+4*norm_vectorsZ(no_points)*ref_track(3,no_points)-2*norm_vectorsX(no_points)*ref_track(1,1)-2*norm_vectorsY(no_points)*ref_track(2,1)-2*norm_vectorsZ(no_points)*ref_track(3,1);

%Solve problem
dev_max_right = ref_track(4,:);
dev_max_left = ref_track(5,:);
dev_max_right(dev_max_right < 0.001)=0.001;
dev_max_left(dev_max_left < 0.001)=0.001;
G = vertcat(eye(no_points),-eye(no_points));
h = [dev_max_right dev_max_left];

alpha = quadprog(H,-f,-G,h);

%% Turning output alpha back into a plotted line
alpha(end+1) = alpha(1);
rate = (alpha + (1/2*(track_width))) / track_width;
rate(rate > 1) = 1;
rate(rate < 0) = 0;
for i = 1:length(ref_track)+1
min_distX(i) = rate(i)*OuterXYZ(i,1) + (1-rate(i))*InnerXYZ(i,1);
min_distY(i) = rate(i)*OuterXYZ(i,2) + (1-rate(i))*InnerXYZ(i,2);
min_distZ(i) = rate(i)*OuterXYZ(i,3) + (1-rate(i))*InnerXYZ(i,3);
end

%Plotting results
% figure
% hold on
% plot3(min_distX, min_distY, min_distZ)
% plot3(InnerX, InnerY, InnerZ)
% plot3(OuterX, OuterY, OuterZ)
% grid on
% axis equal

Min_DistXYZ = [min_distX ; min_distY ; min_distZ];
cubic_spline_inner = cscvn(InnerXYZ(:,[1:end 1])');
cubic_spline_outer = cscvn(OuterXYZ(:,[1:end 1])');
cubic_spline_min_dist = cscvn(Min_DistXYZ(:,[1:end 1]));
figure
hold on
fnplt(cubic_spline_inner,'k',0.5)
fnplt(cubic_spline_outer,'k',0.5)
fnplt(cubic_spline_min_dist,'bl')