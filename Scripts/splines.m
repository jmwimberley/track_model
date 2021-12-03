%%My own one-to-one mapping

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

%%

%% Cubic spline interpolation appending inner and outer
cubic_spline_inner = cscvn(xyzI(:,[1:end 1])');
cubic_spline_outer = cscvn(xyzO(:,[1:end 1])');

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