function plt = plot_corner_pdf(X,varargin)
% plot_corner_pdf.m
% Benjamin Hanson, 2025
% 
%   Given the N-by-D coordinate matrix X and the N-by-1 probability vector 
%   P, creates a D-by-D corner figure. Each grid of the corner PDF contains 
%   a 2D contour plot with isovalues taken to be the 68-95-99.7% isocontours.
% 
% Inputs:
%          X -- N-by-D coordinate matrix
%   varargin -- optional arguments
%               * P -- weights of state vectors (optional)
%               * p -- plotting parameters (optional)
%                   *   color -- isosurface color
%                   * display -- handle visibility
%                   *    name -- display name, if display==1
%                   *   means -- plot weighted mean of point mass PDF
%                   *     axh -- figure axis
%                   *    type -- distribution type
%                   *    hist -- plot histogram on i==j plots
%               * aR -- alpha radius (for 2D and 3D functionality only)
%               * isovalue -- isosurface value(s) to plot (for 2D and 3D functionality only)
%               * axc -- figure axis to copy (because plotting alphaShape
%               automatically changes the aspect ratio and limits of a plot
%               * lbls -- labels of axes
% 
% Example:
% 
%   figure; hold on; 
%   d = 5; mu = randn(1,d); P = randn(d,d); P = P' * P; N = 20; x = zeros(d, N); 
%   for i = 1:d, x(i,:) = linspace(-3 * sqrt(P(i,i)) + mu(i),3 * sqrt(P(i,i)) + mu(i),N); end
%   [X1,X2,X3,X4,X5] = ndgrid(x(1,:)',x(2,:)',x(3,:)',x(4,:)',x(5,:)'); 
%   X = [X1(:) X2(:) X3(:) X4(:) X5(:)]; 
%   lbl = {"$x_1$", "$x_2$", "$x_3$", "$x_4$", "$x_5$"}; p.color = "red"; p.type = "grid"; 
%   P = mvnpdf(X, mu, P); plot_corner_pdf(X,'P',P,'lbls',lbl,'p',p);

% variable arguments - defaults
p.color = {"r", "r", "r"}; 
p.display = 0; 
p.means = 0;
p.axh = gca; 
p.type = "scatter";
aR = [2 2 2]; 
isovalue = [0.68, 0.95, 0.997]; 
axc = p.axh; 

for i=1:2:length(varargin)
    if strcmp('P',varargin{i})
        P = varargin{i+1};
    elseif strcmp('p',varargin{i})
        p = varargin{i+1};
    elseif strcmp('aR',varargin{i})
        aR = varargin{i+1};
    elseif strcmp('isovalue',varargin{i})
        isovalue = varargin{i+1};
    elseif strcmp('axc',varargin{i})
        axc = varargin{i+1}; 
    elseif strcmp('lbls',varargin{i})
        lbls = varargin{i+1}; 
    else
        error(append("Unspecified argument: ", varargin{i}));
    end
end

if ~isfield(p, 'type')
    p.type = "scatter"; 
end

[N, D] = size(X); 
if strcmp(p.type, "scatter")
    P = zeros(N,1);
end
[NP, ~] = size(P); 

% checks and balances
if N~=NP
    error("Incongruous state vector/weight sets.")
end
if ~isfield(p,'color')
    for i = 1:numel(isovalue)
        p.color{i}="r";
    end
else
    if (isstring(p.color))||(ischar(p.color))||((all(size(p.color) == [1,3]))&&(~iscell(p.color)))
        col = p.color; p.color = {}; 
        for i = 1:numel(isovalue)
            p.color{i}=col;
        end
    end 
end
if ~isfield(p,'display')
    p.display=0;
end
if(p.display == 1)
    if ~isfield(p,'name')
        p.display=0;
    end
end
if isfield(p,'name')
    if ~isfield(p,'display')
        p.display=1;
    end
end
if ~isfield(p,'means')
    p.means=0;
end
if ~isfield(p,'axh')
    p.axh=gca;
end
if isfield(p,'type')
    if strcmp(p.type, "curve")
    elseif strcmp(p.type, "grid")
    elseif strcmp(p.type, "scatter")
        if ~isfield(p,'marker')
            p.marker='o';
        end
        if ~isfield(p,'ms')
            p.ms=10;
        end
    elseif strcmp(p.type, "ukf")
        mu = zeros(D,1); 
        for j = 1:N
            mu = mu + P(j, 1) .* X(j, :)'; % Mean
        end
        S = zeros(D,D); 
        for j = 1:N
            S = S + P(j, 2) .* ((X(j, :)' - mu)*(X(j, :)' - mu)'); % Covariance
        end
    else
        error("Unsupported type.")
    end
end
if exist("aR", "var")
    if isscalar(aR)
        if aR < 0
            error("alphaRadius may not be negative."); 
        else
            aR = aR.*ones(1,numel(isovalue));
        end
    else
        if numel(aR) ~= numel(isovalue)
            error("Incongruous alphaRadius/isovalue sets.");
        else
            if any(aR < 0)
                error("alphaRadius may not be negative."); 
            end
        end
    end
end
if (max(isovalue)>1)||(min(isovalue)<0)
    error("Isovalue is outside of probability bounds [0,1].")
end
if ~isa(axc, 'matlab.graphics.axis.Axes')
    error("Copy axis must be an axis variable.")
end
if ~exist("lbls", "var")
    for i = 1:D
        lbls{i} = "$x_" + num2str(i) + "$"; 
    end
else
    [~, Dl] = size(lbls);
    if Dl ~= D
        error("Incompatible label array size.");
    end
end

lbls_bool = 0;
tl = findall(gcf,'Type','tiledlayout');
if isempty(tl)
    t = tiledlayout(D,D); % create only if none exists
    t.TileSpacing = 'compact';
    lbls_bool = 1; 
else
    % verify it matches D x D
    if ~(tl.GridSize(1)==D && tl.GridSize(2)==D)
        error("Existing tiledlayout has incompatible dimensions.");
    end
end

count = 1; plt = {}; 
for i = 1:D
    for j = 1:i
        pause(0.2); nexttile((i - 1) * D + j); p.axh = gca; 
        if lbls_bool
            hold on; box on; set(p.axh, 'FontName', 'times', 'FontSize', 14, "LineWidth", 2);
        end
        if i==j
            if p.hist
                plt{count} = histogram(X(:,i), 'Normalization', 'probability', 'FaceColor', p.color{1}, "EdgeColor", "none");
                if (i == 1)
                    ylabel("Probability", 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
                end
            end
        else
            if strcmp(p.type, "scatter")
                plt{count} = scatter(X(:,j), X(:,i), p.ms, p.color{1}, 'filled', p.marker);
            elseif strcmp(p.type, "ukf")
                p.type = "line";
                [plt{count}, shp] = plot_gaussian_ellipsoid(mu([j,i], 1), S([j,i], [j,i]), 'p', p); 
                p.type = "ukf";
            else
                plt{count} = plot_nongaussian_surface(X(:,[j,i]),P,'p',p); 
            end

            pad_frac = 0.1;     % 10% padding
            eps_min = 1e-9;     % minimum half-range to avoid degenerate limits

            if ~lbls_bool
                if strcmp(p.type, "ukf")
                    shp = shp{1};
                    x_all = [p.axh.XLim, shp(1, :)];
                    y_all = [p.axh.YLim, shp(2, :)];
                else
                    x_all = [p.axh.XLim, X(:,j)'];
                    y_all = [p.axh.YLim, X(:,i)'];
                end
            else
                if strcmp(p.type, "ukf")
                    shp = shp{1};
                    x_all = shp(1, :);
                    y_all = shp(2, :);
                else
                    x_all = X(:,j);
                    y_all = X(:,i);
                end
            end

            % Compute symmetric limits for x
            x_min = min(x_all);
            x_max = max(x_all);
            x_center = 0.5 * (x_max + x_min);
            x_half_range = 0.5 * (x_max - x_min);
            % guard against zero range and add padding
            x_half_range = max(x_half_range, eps_min);
            x_half_range = x_half_range * (1 + pad_frac);
            xlim(p.axh, [x_center - x_half_range, x_center + x_half_range]);

            % Compute symmetric limits for y
            y_min = min(y_all);
            y_max = max(y_all);
            y_center = 0.5 * (y_max + y_min);
            y_half_range = 0.5 * (y_max - y_min);
            % guard against zero range and add padding
            y_half_range = max(y_half_range, eps_min);
            y_half_range = y_half_range * (1 + pad_frac);
            ylim(p.axh, [y_center - y_half_range, y_center + y_half_range]);

            % lock the limits so MATLAB won't autoscale them later
            p.axh.XLimMode = 'manual';
            p.axh.YLimMode = 'manual';
        end
        
        if lbls_bool
            if j == 1 && (i ~= 1)
                ylabel(lbls{i}, 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
            end
            if i==D
                xlabel(lbls{j}, 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
            end
            if i ~= D, xticks([]); end
            if (j ~= 1), yticks([]); end
        end
        axis normal; 
        count = count + 1; 
        drawnow; 
    end
end
end