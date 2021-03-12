function PlotPCPhantom(phantom_cell)

K = length(phantom_cell);
n = 100;
m = 128;

% Create colorscale blue to red
cond_data = zeros(K,1);
for i = 1:K
    cond_data(i)= phantom_cell{i}.amplitude;
end
maxC = max(cond_data);
minC = min(cond_data);
C = redblue(n);
a = (maxC-minC)/(n-1);
b = minC-a;

tight_subplot(1,1,[0,0],[-0.2,-0.2],[-0.2,-0.2]);

for i = 1:K
    hold on;
    if isa(phantom_cell{i},'PcEllipsoidConductivity')
        if a ~= 0
            c = ceil((phantom_cell{i}.amplitude-b)/a);
        else
            c = n;
        end
        [X,Y,Z]=ellipsoid(0,0,0,phantom_cell{i}.radii(1),phantom_cell{i}.radii(2),phantom_cell{i}.radii(3),m);
        surf(phantom_cell{i}.center(1)+phantom_cell{i}.axes(1,1)*X+phantom_cell{i}.axes(2,1)*Y+phantom_cell{i}.axes(3,1)*Z, ...
            phantom_cell{i}.center(2)+phantom_cell{i}.axes(1,2)*X+phantom_cell{i}.axes(2,2)*Y+phantom_cell{i}.axes(3,2)*Z, ... ,
            phantom_cell{i}.center(3)+phantom_cell{i}.axes(1,3)*X+phantom_cell{i}.axes(2,3)*Y+phantom_cell{i}.axes(3,3)*Z, ...
            'EdgeColor','none','FaceColor',C(c,:));
        
    elseif isa(phantom_cell{i},'PcBallConductivity')
        if a ~= 0
            c = ceil((phantom_cell{i}.amplitude-b)/a);
        else
            c = n;
        end
        [X,Y,Z]=sphere(m);
        surf(phantom_cell{i}.center(1)+phantom_cell{i}.radius*X,phantom_cell{i}.center(2)+phantom_cell{i}.radius*Y,phantom_cell{i}.center(3)+phantom_cell{i}.radius*Z,'EdgeColor','none','FaceColor',C(c,:));
    else
        error('Please construct conductivity data from either ellipsoids or spheres');
    end
end


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