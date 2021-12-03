%% LOAD DATA
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

%%CONVERSION TO X Y COORDS
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

%% METHOD 1 - ARC LENGTH AVERAGING
xyzInner = [xI yI RelInnerAlt];
xyzOuter = [xO yO RelOuterAlt];

xyzI = interparc(10000,xI,yI,RelInnerAlt);
xyzO = interparc(10000,xO,yO,RelOuterAlt);
xyzC = (xyzI+xyzO)/2;

%%Plotting the 3 interparc splines
hold on
plot3(xyzI(:,1),xyzI(:,2),xyzI(:,3))
plot3(xyzO(:,1),xyzO(:,2),xyzO(:,3))
plot3(xyzC(:,1),xyzC(:,2),xyzC(:,3))
grid on
axis equal

%% METHOD 2 - COEFFICIENTS OF CUBIC SPLINE AVERAGING

for j=2:length(xyzOuter)
    for i=2:length(xyzInner)
        xOuter(1,1) = xO(end-1);
        yOuter(1,1) = yO(end-1);
        zOuter(1,1) = RelOuterAlt(end-1);
        xOuter(i,1) = (((i-1)/(length(xyzInner)))*(xO(1)-xO(end-1)))+xO(end-1);
        yOuter(i,1) = (((i-1)/(length(xyzInner)))*(yO(1)-yO(end-1)))+yO(end-1);
        zOuter(i,1) = (((i-1)/(length(xyzInner)))*(RelOuterAlt(1)-RelOuterAlt(end-1)))+RelOuterAlt(end-1);
        xOuter(1,j) = xO(j-1);
        yOuter(1,j) = yO(j-1);
        zOuter(1,j) = RelOuterAlt(j-1);
        xOuter(i,j) = (((i-1)/(length(xyzInner)))*(xO(j)-xO(j-1)))+xO(j-1);
        yOuter(i,j) = (((i-1)/(length(xyzInner)))*(yO(j)-yO(j-1)))+yO(j-1);
        zOuter(i,j) = (((i-1)/(length(xyzInner)))*(RelOuterAlt(j)-RelOuterAlt(j-1)))+RelOuterAlt(j-1);
    end
end

for j=2:length(xyzInner)
    for i=2:length(xyzOuter)
        xInner(1,1) = xI(end-1);
        yInner(1,1) = yI(end-1);
        zInner(1,1) = RelInnerAlt(end-1);
        xInner(i,1) = (((i-1)/(length(xyzOuter)))*(xI(1)-xI(end-1)))+xI(end-1);
        yInner(i,1) = (((i-1)/(length(xyzOuter)))*(yI(1)-yI(end-1)))+yI(end-1);
        zInner(i,1) = (((i-1)/(length(xyzOuter)))*(RelInnerAlt(1)-RelInnerAlt(end-1)))+RelInnerAlt(end-1);
        xInner(1,j) = xI(j-1);
        yInner(1,j) = yI(j-1);
        zInner(1,j) = RelInnerAlt(j-1);
        xInner(i,j) = (((i-1)/(length(xyzOuter)))*(xI(j)-xI(j-1)))+xI(j-1);
        yInner(i,j) = (((i-1)/(length(xyzOuter)))*(yI(j)-yI(j-1)))+yI(j-1);
        zInner(i,j) = (((i-1)/(length(xyzOuter)))*(RelInnerAlt(j)-RelInnerAlt(j-1)))+RelInnerAlt(j-1);
    end
end
OuterEqual = [xOuter(:) yOuter(:) zOuter(:)];
InnerEqual = [xInner(:) yInner(:) zInner(:)];

%%Cubic spline interpolation appending inner and outer
cubic_spline_inner = cscvn(InnerXYZ(:,[1:end 1])');
cubic_spline_outer = cscvn(OuterXYZ(:,[1:end 1])');

%Building identical centre structure from inner and outer data
cubic_spline_centre.form = 'pp';
cubic_spline_centre.breaks = (cubic_spline_inner.breaks+cubic_spline_outer.breaks)/2;
cubic_spline_centre.coefs = (cubic_spline_inner.coefs+cubic_spline_outer.coefs)/2;
cubic_spline_centre.pieces = length(cubic_spline_centre.breaks)-1;
cubic_spline_centre.order = 4;
cubic_spline_centre.dim = 4;

hold on
fnplt(cubic_spline_inner,'k',1)
fnplt(cubic_spline_outer,'k',1)
fnplt(cubic_spline_centre,'r',1)
%Plot those bad boys

%% METHOD 3 - INDEXING

for k = 1 : length(OuterXYZ)
    distances = sqrt((xO(k) - xI) .^ 2 + (yO(k) - yI) .^ 2 + (RelOuterAlt(k) - RelInnerAlt) .^ 2);
    [minDistance, indexOfMin] = min(distances);
    xCentre(k) = mean([xO(k), xI(indexOfMin)]);
    yCentre(k) = mean([yO(k), yI(indexOfMin)]);
    zCentre(k) = mean([RelOuterAlt(k), RelInnerAlt(indexOfMin)]);
end

%% PLOTTING ALL 3

hold on
fnplt(cubic_spline_inner,'k',0.5)
fnplt(cubic_spline_outer,'k',0.5)
fnplt(cubic_spline_centre,'g')
plot3(xyzC(:,1),xyzC(:,2),xyzC(:,3))
plot3(xCentre,yCentre,zCentre,'g')
axis equal