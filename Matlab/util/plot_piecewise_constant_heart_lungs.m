function plot_piecewise_constant_heart_lungs(phantomsize)
% PLOT_PIECEWISE_CONSTANT_HEART_LUNGS(size)
% Plots a piecewise conductivity phantom whose 2D view looks like a heart
% and lungs.
%
% phantomsize is either 'small' or 'big'.
% If not specified, the size is put to 'small'.

% Phantom size
if nargin==0 % small size
    k=3/5;
elseif nargin==1
    if strcmp(phantomsize,'small')
        k=3/5;
    elseif strcmp(phantomsize,'big')
        k=1;
    else
        error('Bad input argument, choose nothing, small or big');
    end
else
    error('Bad number of input arguments');
end

% Ratio between radii of the prolate spheroids
h=2;
k =1;

%------- Parameters ---------% 
% inner extent
rotang=5*pi/12;
% factor in radius
k1 = 0.65;
k2 = 1.3;
% Prolate spheroid on the right
% Center
xr0=sin(rotang)*0.45;
yr0=cos(rotang)*0.45;
zr0=0;
% Radius
dr=k1*0.42;
% Conductivity
ar=0.5;

% Prolate spheroid on the left
xl0=-sin(rotang)*0.55;
yl0=cos(rotang)*0.55;
zl0=0;
% Smallest radius (up to the product by k)
dl=k1*0.36;
% Conductivity
al=0.5;

% Ball
% Center
xb0=-0.09;
yb0=-0.55;
zb0=0;
% Radius
db=k2*0.21;
% Conductivity
ab=2;

tight_subplot(1,1,[0,0],[-0.2,-0.2],[-0.2,-0.2]);
hold on;
[X,Y,Z]=ellipsoid(0,0,0,k*h*dr,k*dr,k*dr,128);
hr=surf(xr0+cos(rotang)*X+sin(rotang)*Y,yr0-sin(rotang)*X+cos(rotang)*Y,zr0+Z,'EdgeColor','none','FaceColor','blue');
[X,Y,Z]=ellipsoid(0,0,0,k*h*dl,k*dl,k*dl,128);
hl=surf(xl0+cos(rotang)*X-sin(rotang)*Y,yl0+sin(rotang)*X+cos(rotang)*Y,zl0+Z,'EdgeColor','none','FaceColor','blue');
[X,Y,Z]=sphere(128);
hb=surf(xb0+k*db*X,yb0+k*db*Y,zb0+k*db*Z,'EdgeColor','none','FaceColor','red');

title('');
axis([-1,1,-1,1,-1,1]);axis equal;
view(3);
grid off;
camlight headlight;
[X,Y,Z]=sphere(24);
h=mesh(X,Y,Z);
set(h,'facecolor','none','edgecolor','k','linestyle',':');
hold off;
set(gcf,'color',[1 1 1]);
axis off;