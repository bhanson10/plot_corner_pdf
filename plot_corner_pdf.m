function plt = plot_corner_pdf(X,P,lbl,Xtruth,isovalue,p)
%Plot Corner PDF
%   PLT = PLOT_CORNER_PDF(X,P) takes the N-by-D coordinate matrix X and the
%   D-by-1 probability vector P and creates a D-by-D corner figure. Each grid
%   of the corner PDF contains a 2D contour plot with isovalues taken to be
%   the 68-95-99.7% isocontours. 
%
%   PLT = PLOT_CORNER_PDF(X,P,LBL) takes the N-by-D coordinate matrix X 
%   and the D-by-1 probability vector P and creates a D-by-D corner figure.
%   Each grid of the corner PDF contains a 2D contour plot with isovalues 
%   taken to be the 68-95-99.7% isocontours. The 1-by-D cell array LBL adds
%   the axes lables to the left y-axes and the bottom x-axes.
%
%   PLT = PLOT_CORNER_PDF(X,P,LBL,XTRUTH) takes the N-by-D coordinate matrix 
%   X  and the D-by-1 probability vector P and creates a D-by-D corner figure.
%   Each grid of the corner PDF contains a 2D contour plot with isovalues 
%   taken to be the 68-95-99.7% isocontours. The 1-by-D cell array LBL adds
%   the axes lables to the left y-axes and the bottom x-axes. The D-by-1 truth 
%   vector XTRUTH plots the true coordinate on each contour plot; on (i,j) 
%   plots when i=j, this is represented as  a vline, and on (i,j) when i~=j
%   this is represented as a scatter.
%
%   PLT = PLOT_CORNER_PDF(X,P,LBL,XTRUTH,ISOVALUE) takes the N-by-D coordinate 
%   matrix X and the D-by-1 probability vector P and creates a D-by-D corner
%   figure. Each grid of the corner PDF contains a 2D contour plot with 
%   1-by-M ISOVALUE isocontours. The 1-by-D cell array LBL adds the axes 
%   lables to the left y-axes and the bottom x-axes. The D-by-1 truth vector 
%   XTRUTH plots the true coordinate on each contour plot. 
%
%   PLT = PLOT_CORNER_PDF(X,P,LBL,XTRUTH,ISOVALUE,p) takes the N-by-D 
%   coordinate matrix X and the D-by-1 probability vector P and creates a 
%   D-by-D corner figure. Each grid of the corner PDF contains a 2D contour 
%   plot with 1-by-M ISOVALUE isocontours. The 1-by-D cell array LBL adds 
%   the axes lables to the left y-axes and the bottom x-axes. The D-by-1 
%   truth vector XTRUTH plots the true coordinate on each contour plot. 
%   Each grid of the corner PDF contains a 2D contour plot with ISOVALUE 
%   isocontours. p is used for plotting parameters. 
%       * p -- plotting parameters (optional)
%            *   color -- isosurface color
%            * display -- handle visibility
%            *    name -- display name, if display==1
%            *   means -- plot weighted mean of point mass PDF
%            *     axh -- figure axis
%            *   alpha -- surface visibility
%            *    type -- distribution type
%
%   Requirements: 
%      - X must be cartesian gridded coordinates
%      - plot_nongaussian_surface.m
%
%   Example: 
%     figure; hold all; 
%     d = 5; mu = randn(1,d); P = randn(d,d); P = P' * P; N = 20; x = zeros(d, N); 
%     for i = 1:d, x(i,:) = linspace(-3 * sqrt(P(i,i)) + mu(i),3 * sqrt(P(i,i)) + mu(i),N); end
%     [X1,X2,X3,X4,X5] = ndgrid(x(1,:)',x(2,:)',x(3,:)',x(4,:)',x(5,:)'); 
%     X = [X1(:) X2(:) X3(:) X4(:) X5(:)]; 
%     lbl = {"$x_1$", "$x_2$", "$x_3$", "$x_4$", "$x_5$"}; p.color = "jet"; 
%     P = mvnpdf(X, mu, P); plot_corner_pdf(X,P,lbl,[],[], p);
% 
% Copyright 2025 by Benjamin L. Hanson, published under BSD 2-Clause License. 

if nargin==2
    [n, d] = size(X); 
    if d < 2
        error("BadXDimension");
    end
    [n_P, ~] = size(P);
    if n_P ~= n
        error("IncompatiblePSize");
    end
    lbl = {}; 
    Xtruth = [];
    isovalue = [normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)];
elseif nargin==3
    [n, d] = size(X); 
    if d < 2
        error("BadXDimension");
    end
    [n_P, ~] = size(P);
    if n_P ~= n
        error("IncompatiblePSize");
    end
    [~, d_lbl] = size(lbl); 
    if d_lbl ~= d
        error("IncompatibleLblSize");
    end
    Xtruth = [];
    isovalue = [normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)];
elseif nargin==4
    [n, d] = size(X); 
    if d < 2
        error("BadXDimension");
    end
    [n_P, ~] = size(P);
    if n_P ~= n
        error("IncompatiblePSize");
    end
    [~, d_lbl] = size(lbl); 
    if d_lbl ~= d
        error("IncompatibleLblSize");
    end
    if ~isempty(Xtruth)
        [d_t, ~] = size(Xtruth); 
        if d_t ~= d
            error("IncompatibleXtruthSize");
        end
    end
    isovalue = [normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)];
elseif nargin==5||nargin==6
    [n, d] = size(X); 
    if d < 2
        error("BadXDimension");
    end
    [n_P, ~] = size(P);
    if n_P ~= n
        error("IncompatiblePSize");
    end
    [~, d_lbl] = size(lbl); 
    if d_lbl ~= d
        error("IncompatibleLblSize");
    end
    if ~isempty(Xtruth)
        [d_t, ~] = size(Xtruth); 
        if d_t ~= d
            error("IncompatibleXtruthSize");
        end
    end
    if ~isempty(isovalue)
        [~, M] = size(isovalue);
        if M < 3
            error("MoreIsovaluesRequired");
        end
    else
        isovalue = [normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)];
    end
else
    error("BadNumberOfInputs");
end

tiledlayout(d,d);

count = 1; plt = {}; 
for i = 1:d
    for j = 1:i
        nexttile((i-1)*d + j); hold on;
        if i==j
            if exist('p','var')
                plt{count} = plot_nongaussian_surface(X(:,i),P,isovalue,p); 
            else
                plt{count} = plot_nongaussian_surface(X(:,i),P,isovalue); 
            end
            if ~isempty(Xtruth)
                xline(Xtruth(i), "r-", "LineWidth", 2);
            end
        else
            if exist('p','var')
                plt{count} = plot_nongaussian_surface(X(:,[j,i]),P,isovalue,p); 
            else
                plt{count} = plot_nongaussian_surface(X(:,[j,i]),P,isovalue); 
            end
            if ~isempty(Xtruth)
                scatter(Xtruth(i), Xtruth(j), 50, "r", "filled", "Marker", "d");
            end
        end

        if j == 1
            ylabel(lbl{i}, 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
        elseif i ~= j
            yticklabels([]);
        end
        if i==d
            xlabel(lbl{j}, 'FontSize', 18, 'FontName', 'Times', 'Interpreter', 'latex');
        else
            xticklabels([]);
        end
    end
end
end