function map = haxby_hk(m)

%Get default colormap size if user doesn't specify one
if nargin < 1
    m = size(get(gcf,'colormap'),1);
end

%Haxby color palette from GMT
c = [10  0 121;
     0   10  200;
     26  102 240;
     50  190 255;
     106 235 225;
     205 255 162;
     247 215 104;
     244 117 75 ;
     255 124 124;
     255 196 196;
     255 254 253;];

%Set a color if m = 1;
if m==1
    map = [c(1,1) c(1,2) c(1,3)]/255;
%Interpolate if m > 1;
else
    ncolors = 11;
    p = 1:(m-1)/(ncolors-1):m;
    r = interp1(p,c(:,1),1:m);
    g = interp1(p,c(:,2),1:m);
    b = interp1(p,c(:,3),1:m);
    map = [r' g' b']/255;
end
