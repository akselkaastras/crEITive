function [ha, pos] = Plot3Dc(reconid,dnmapid,xslice,yslice,zslice,h,span, varargin)
% Mandatory input:
%   type    = 'conductivity'
%   reconid - number or vector of ids. If zero we must have dnmapid ? 0
%   dnmapid - dnmapid of true conductivity which user would like to be
%             plotted. If dnmapid = 0, then no true conductivity is plotted
%   xslice  - vector or scalar representing one or more slices orthogonal
%             to x-plane
%   yslice  - vector or scalar representing one or more slices orthogonal
%             to y-plane
%   zslice  - vector or scalar representing one or more slices orthogonal
%             to z-plane
%   h       - fineness of grid on which type is interpolated
%   span    - [a,b] determines the span of the colorbar
%   varargin- If type is complex one can choose real or imagainary parts by 
%             'real' or 'imag' resp. (default is real)
%
%
% Output:
%   ha      - array of handles of the axes objects starting from upper 
%             left corner, going row-wise as in subplot
%   pos     - positions of the axes objects

type = 'conductivity';
if dnmapid == 0 && sum(reconid) == 0
    error('Input something to plot!');
end

% check varargin and how to plot complex
cmplx = 'real';
if nargin == 8
    cmplx = varargin{1};
end
% initialize index offset if true conductivity is to be plotted
if dnmapid == 0
    k = 0;
else
    k = 1; % true conductivity is to be plotted
end

M = length(reconid);
if sum(reconid) == 0
    M = M-1;
end

N = k+M;
for i = M:-1:1
    if ~reconid
        break;
    end
    
    [filename,info] = getfilename_reconid(reconid(i),type);

    % read file depending on type
    if strcmp(type,'conductivity')
        [x,y,z,c] = read_conductivity(filename);
    elseif strcmp(type,'qhat')
        [x,y,z,c] = read_qhat(filename,'meshgrid');
    elseif strcmp(type,'q')
        [x,y,z,c] = read_q(filename);
    else
        error('Please choose a type');
    end

    % insert to struct:
    cond(i+k).radial = 0;
    cond(i+k).true = 0;
    % Get parameter values
    cond(i+k).reconid = info{2};
    cond(i+k).reconmethod = info{3};
    cond(i+k).nd = info{5};
    cond(i+k).zeta = info{6};
    cond(i+k).ift = info{7};
    cond(i+k).ngrid = info{9};
    cond(i+k).pkappa = info{11};
    dnmapdat = split(info{13},'.');
    cond(i+k).dnmapid = dnmapdat{1};
    cond(i+k).x = x;
    cond(i+k).y = y;
    cond(i+k).z = z;
    cond(i+k).c = c;
end

% Deal with true conductivities from dnmap here:
if dnmapid ~=0
    [filename,radial] = getfilename_dnmapid(dnmapid,type);

    if radial % real radially symmetric conductivity
        if strcmp(type,'conductivity')
            [r,c] = read_conductivity(filename);
        elseif strcmp(type,'qhat')
            [r,c] = read_qhat(filename,'meshgrid');
        elseif strcmp(type,'q')
            [r,c] = read_q(filename);
        else
            error('Please choose a type');
        end

        % Radial mesh
        phi = 0:h:(2*pi+h);
        [R, PHI] = meshgrid(r,phi);
        % interpolate
        F = interp1(r,c,R);

        % insert to struct:
        cond(1).x = R.*cos(PHI);
        cond(1).y = R.*sin(PHI);
        cond(1).radial = 1;
        cond(1).c = F;
        cond(1).true = 1;

    else
        if strcmp(type,'conductivity')
            [x,y,z,c] = read_conductivity(filename);
        elseif strcmp(type,'qhat')
            [x,y,z,c] = read_qhat(filename,'meshgrid');
        elseif strcmp(type,'q')
            [x,y,z,c] = read_q(filename);
        else
            error('Please choose a type');
        end

        % insert to struct:
        cond(1).x = x;
        cond(1).y = y;
        cond(1).z = z;
        cond(1).c = c;
        cond(1).radial = 0;
        cond(1).true = 1;
    end
end


% Determine arrangement of plots and width and height
if N == 1
    Nh = 1; % horizontal
    Nw =1; % vertical
    width = 12;
    height = 11;
    gap = [0, 0];
    margh = [0.06, 0.05];
    margw = [0.03, 0];
elseif rem(N,2) == 0
    Nh = N/2; % horizontal
    Nw = 2; % vertical
    width = 22;
    height = ceil(N/2)*10;
    gap = [0.1, 0.05];
    margh = [0.05, 0.05];
    margw = [0.05, 0.15];
else
    error('Please choose {1,2,4,6,...} plots')
end


% Tight subplot
%figure;
[ha,pos] = tight_subplot(Nh, Nw,gap,margh,margw);

% Interpolate and plot
for ii = 1:N
    axes(ha(ii));
    
    % deal with radial conductivities here
    if cond(ii).radial
        surf(cond(ii).x,cond(ii).y,cond(ii).c,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
        titletext = 'Conductivity, \gamma';
        ax = gca;
        figsettings(titletext,span,ax)
    else %non radial here
        if strcmp(type,'conductivity')
            k1 = length(xslice);
            k2 = length(yslice);
            k3 = length(zslice);
            SLICE = [xslice yslice zslice];
            
            for o = (k1+k2+k3):-1:1
                if o <= k1
                    dimension = 1;
                elseif o > k1 && o<= k2+k1
                    dimension = 2;
                else
                    dimension = 3;
                end
                [xs,ys,zs,P,Q] = cross_plane(-1,1,h,SLICE(o),dimension);
                % interpolate
                Vq = scatteredInterpolant(cond(ii).x,cond(ii).y,cond(ii).z,cond(ii).c,'linear','none');
                Vq = Vq(xs,ys,zs);
                
                %thress = (Vq < 0.5 | Vq > 2);
                %thress = xs.^2 + ys.^2 + zs.^2 > 1;
                %Vq(thress) = NaN;

                % if complex
                if ~isreal(Vq)
                    disp('Plot is complex...');
                    if strcmp(cmplx,'1norm')
                        Vq = abs(real(Vq))+abs(imag(Vq));
                        disp('Plotting 1-norm');
                    elseif strcmp(cmplx,'real')
                        Vq = real(Vq);
                        disp('Plotting the real part');
                    elseif strcmp(cmplx,'imag')
                        Vq = imag(Vq);
                        disp('Plotting the imaginary part');
                    end
                end

                % plot
                surf(xs,ys,zs,Vq);
                hold on
            end
            hold off
            
            titletext = 'Conductivity, $$\gamma$$';

            ax = gca;
            figsettings(titletext,span,[-1,1],ax);
            
            % deal with very high values
            if max(max(Vq))>2*span(2)
                caxis auto;
            end
            
        end  
    end
end

% Determine fitting width and height of figure
set(gcf,'units','centimeters','position',[10,10,width,height])

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Additional functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function figsettings(titletext,span,xyspan,ax)
    view(3);
    shading('interp');
    axis image;
    caxis manual;
    caxis(span)
    colormap jet;
    axis([-1,1,-1,1,-1,1]);axis equal;
    view([-30,25]);
    %axis off;
    %colorbar;
    grid off
    a = xyspan(1);
    b = xyspan(2);
    xyspan = linspace(ceil(a),floor(b),5);
    xlim([a,b])
    ylim([a,b])
    yticks(xyspan)
    xticks(xyspan)
    ax.FontSize = 12; 
    box on;
    ax.LineWidth = 1.5;
    title(titletext, 'Interpreter', 'Latex','FontSize',14)
    %set(gca,'visible','off')
    %set(gca,'xtick',[])
end
function [filename,info] = getfilename_reconid(reconidid,type)
filenameshort = sprintf('results/%s/%s_%d_*',type,type,reconidid);
matches = struct2cell(dir(filenameshort));
s2 = size(matches,2);
% check if reconid is unique
if s2 > 1
    error('Nonunique reconstruction id, = %d. Please choose a unique id', reconidid);
elseif s2 == 0
    error('No such reconstruction id, reconid = %d, exists', reconidid);
end
% get info
info = regexp(matches{1},'_','split');
filename = sprintf('results/%s/%s',type,matches{1});
end
function [filename,radial] = getfilename_dnmapid(dnmapidid,type)

filenameshort = sprintf('results/%s/%s_dnmapid_%d_*',type,type,dnmapidid);
matches = struct2cell(dir(filenameshort));
filename = sprintf('results/%s/%s',type,matches{1});
s2 = size(matches,2);
% check if reconid is unique
if s2 > 1
    error('Nonunique dnmap id, = %d. Please choose a unique id', dnmapidid);
elseif s2 == 0
    error('No such dnmap id, = %d, exists', dnmapidid);
end
% check if radial
fid = fopen(filename,'r');

% Check error
if fid == -1
    error(['Error when opening file ',filename]);
end

% Read the first line and remove leading and trailing white space
firstline=strtrim(fgetl(fid));
radial = strcmp(firstline,'radial');
end
function [xs,ys,zs,P,Q] = cross_plane(a,b,h,x0,dimension)
% Crossectional plane determined by dimension and x0
if dimension == 1
    v = [1,0,0];
elseif dimension == 2
    v = [0,1,0];
else
    v = [0,0,1];
end
v = v/norm(v);
w = null(v);
[P,Q] = meshgrid(a:h:b);
x1 = zeros(3,1);
x1(dimension) = x0;
xs = w(1,1)*P+w(1,2)*Q+x1(1); % x-coordinates of slice
ys = w(2,1)*P+w(2,2)*Q+x1(2); % y-coordinates of slice
zs = w(3,1)*P+w(3,2)*Q+x1(3); % z-coordinates of slice
end