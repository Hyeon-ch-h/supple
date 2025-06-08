clear; clc;

%% (0) Read data


% Ask user for the input filename (without .csv)
filename = "99 Circles 1 30_grad";
cd = 250;         % diameter [um]
x0 = 2;
y0 = 3.125;      % offset from origin? [mm]
rotationAngle = 30;  % rotation in degree
choice = ['BdB'];     % 'B','BdB','dB2','Sm','BdB_l' / choose one

% Number of Circles and Array Dimensions
numCircles = 4;  % Number of Circles
rows = 2;  % Number of rows
cols = 2;  % Number of columns


%% (1) Seprate Data
% Define full file paths
input_file = filename + ".csv";            % Input CSV file
% Read the CSV file (skipping header)
loaded= readtable(input_file);
X = loaded{:, 1} - x0; 
Y = loaded{:, 2} - y0;  
Bx = loaded{:,3};
By = loaded{:,4};

% Z depending on choice
if strcmp(choice, 'B')
    Z = loaded{:, 5}; 
    xname = '|B| [T]';  % x label on histogram
    color_limit=[0, 1.0];
    ll='linear';
elseif strcmp(choice, 'BdB') || strcmp(choice, 'BdB_l')
    Z = loaded{:, 6}/2; 
xname = '$\mathbf{|\left( B \cdot \nabla \right) B|} \quad \mathbf{[T^2/m]}$';
    color_limit=[1, 10^7];
    ll='log';
elseif strcmp(choice, 'dB2')
    Z = loaded{:, 6}; 
xname = '$\nabla B^2 \quad [T^2/m]$';
    color_limit=[1, 10^7];
    ll='log';
elseif strcmp(choice, 'Sm')
    Z = loaded{:, 6}/2/(4*pi*10^(-7)); 
    xname = 'S_m [TA/m^2]';
    color_limit=[1, 10^7];
    ll='log';
else
    error('Invalid choice. Please select a valid option.');
end

%% (2) Circle center
r= cd/1000/2;

% Circle-to-circle spacing (based on radius only, with no extra margin)
gridSpacing = 2 * r; % Set spacing to 2 × radius so circles are just touching

% Generate circle center coordinates in a grid layout
[xGrid, yGrid] = meshgrid(linspace(-gridSpacing * (cols - 1) / 2, gridSpacing * (cols - 1) / 2, cols), ...
                          linspace(-gridSpacing * (rows - 1) / 2, gridSpacing * (rows - 1) / 2, rows));

% Combine coordinates into a single array (numCircles × 2)
circleCenters = [xGrid(:), yGrid(:)]; 

% Convert rotation angle from degrees to radians
theta = deg2rad(rotationAngle); 

% Apply rotation matrix to all circle centers
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
circleCenters_rotated = (R * circleCenters')';


%% (3) Construct 2D Grid
xUni = unique(X);
yUni = unique(Y);
Nx = length(xUni);
Ny_ = length(yUni);

[Xq, Yq] = meshgrid(xUni, yUni);  % size [Ny_, Nx]

Z2D = nan(Ny_, Nx);  
Bx2D=nan(Ny_, Nx);
By2D=nan(Ny_, Nx);

for k = 1:length(X)
    [~, ix] = min(abs(xUni - X(k)));
    [~, iy] = min(abs(yUni - Y(k)));
    Z2D(iy, ix) = Z(k);
    Bx2D(iy, ix) = Bx(k);
    By2D(iy, ix) = By(k);
end

%% (4) Mask for Excluding Circle Interior/Boundary (2D)
excludeMask2D = false(Ny_, Nx);

for i = 1:numCircles
    cx = circleCenters_rotated(i,1);
    cy = circleCenters_rotated(i,2);

    excludeMask2D = excludeMask2D | ...
        ((Xq - cx).^2 + (Yq - cy).^2 <= r^2);
end
remainMask2D = ~excludeMask2D;


%% Check Circle Centers
%

% (4-1) Create an outer rectangle using the centers of the outermost circles
% 🔹 Find min/max X and Y from the already rotated circle array (based on unrotated state)
minX = min(circleCenters(:,1)); maxX = max(circleCenters(:,1));
minY = min(circleCenters(:,2)); maxY = max(circleCenters(:,2));

% 🔹 Define rectangle corners before rotation (based on unrotated grid)
rectCorners = [
    minX, minY;
    maxX, minY;
    maxX, maxY;
    minX, maxY
];

% 🔹 Rotate the rectangle by the same angle as the circles
rectCorners_rotated = (R * rectCorners')';

% 🔹 Use `inpolygon` to remove points outside the rotated rectangle
outerMask = inpolygon(Xq, Yq, rectCorners_rotated(:,1), rectCorners_rotated(:,2)) == 0;

% 🔹 Extract closed regions by excluding both the circle interior and outside the rectangle
closedRegionMask = ~excludeMask2D & ~outerMask;

% (4-2) Visualization (plot)
%{
figure;
hold on;
axis equal;
grid on;
xlabel('X-axis');
ylabel('Y-axis');
title(sprintf('Rotated vs. Original Circles & Rectangle (Angle = %d°)', rotationAngle));

% 🔵 Plot original circle centers
plot(circleCenters(:,1), circleCenters(:,2), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Original Circles');

% 🔴 Plot rotated circle centers
plot(circleCenters_rotated(:,1), circleCenters_rotated(:,2), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Rotated Circles');

% 🟢 Plot boundary of original rectangle
plot([rectCorners(:,1); rectCorners(1,1)], ...
     [rectCorners(:,2); rectCorners(1,2)], ...
     'g-', 'LineWidth', 1.5, 'DisplayName', 'Original Rectangle');

% 🟣 Plot boundary of rotated rectangle
plot([rectCorners_rotated(:,1); rectCorners_rotated(1,1)], ...
     [rectCorners_rotated(:,2); rectCorners_rotated(1,2)], ...
     'm-', 'LineWidth', 1.5, 'DisplayName', 'Rotated Rectangle');

% 🔹 Mark the rotation center (0, 0)
plot(0, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Rotation Center');

legend;
hold off;

%}

%% (5) Histogram + save (CSV)
% (1) Entire area
entireVals = Z2D(~isnan(Z2D));

% (2) Closed regions formed by touching circles
closedRegionVals = Z2D(closedRegionMask);
closedRegionVals = closedRegionVals(~isnan(closedRegionVals));


% --- Handling non-positive values (for log scale plots) ---
if strcmp(choice, 'B')
    entireValsPos   = entireVals;
    closedRegionValsPos   = closedRegionVals;
else
    entireValsPos   = entireVals(entireVals>0);
    closedRegionValsPos   = closedRegionVals(closedRegionVals>0);
end

%%
% --- Set histogram bins ---
numBins = 101;

if strcmp(choice, 'B') %|| strcmp(choice, 'BdB_l')
    minV=0;
    maxV=2;
    edges = linspace(minV, maxV, numBins);   % Linear binning
    % Arithmetic mean of linear bin edges
    binCenters = (edges(1:end-1) + edges(2:end)) / 2;
else
    minV=10^-6;
    maxV=10^9;
    edges = logspace(log10(minV), log10(maxV), numBins);   % Logarithmic binning
    % Geometric mean of logarithmic bin edges
    binCenters = sqrt(edges(1:end-1) .* edges(2:end));
end

% Start and end values for each bin
binStart = edges(1:end-1);
binEnd   = edges(2:end); 


% Set color range for turbo colormap (fixed between 1 and 10^7)
numBins = numel(binCenters);
colormapData = turbo(numBins);   % Generate turbo colormap

% Normalize bin centers using log scale (adjust values < 1)
binCenters_fixed = max(binCenters, 1);  % Replace values < 1 with 1 for log transform
log_binCenters = log10(binCenters_fixed);  

log_min = log10(1);      % Minimum of colormap (log(1))
log_max = log10(10^7);   % Maximum of colormap (log(10^7))
log_binCenters = (log_binCenters - log_min) / (log_max - log_min);  % Normalize to [0, 1]

% Adjust color indices to prevent exceeding colormap range
colorIndices = round(1 + log_binCenters * (numBins - 1));  
colorIndices = min(colorIndices, numBins);  




%% Plot
figure1 = figure('Name','Histograms: Entire / Excluded / Remain','NumberTitle','off');
set(gcf, 'Units', 'normalized', 'OuterPosition', [0,0,0.6,0.8]);

if strcmp(choice, 'BdB_l')
set(gcf, 'Units', 'normalized', 'OuterPosition', [0,0,0.7,0.5]);
end

tiledlayout(1,2, 'TileSpacing', 'Compact', 'Padding', 'Compact');   % Create 2 tiles

% (a) Entire area
nexttile; % First histogram
hold on;

% Compute histogram data
[counts_entire, binEdges] = histcounts(entireValsPos, edges);

% Compute bin centers and widths
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
binWidths = diff(binEdges);  

cumulative_entire = cumsum(counts_entire) / sum(counts_entire);

for i = 1:numBins
    b = bar(binCenters(i), counts_entire(i), 'FaceColor', colormapData(colorIndices(i), :), ...
        'EdgeColor', 'k', 'BarWidth', binWidths(i));
    if strcmp(choice, 'BdB_l')
        b.LineWidth = 1;
    end
end

% Set X-axis scale
if strcmp(choice, 'B') %|| strcmp(choice, 'BdB_l')
    set(gca, 'XScale','linear');
else
    set(gca, 'XScale','log');
    xlim([2*10^-7, 10^10]);
end


% Add cumulative line on right Y-axis
yyaxis right;
set(gca, 'YColor', 'b');
plot(binCenters, cumulative_entire, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Cumulative');
ylabel('Cumulative Probability','color','b', 'FontWeight', 'bold');


% Set left Y-axis to log scale
yyaxis left;
set(gca, 'YScale','log');
    ylim([1, 10^6]);
if strcmp(choice, 'BdB_l')
    set(gca, 'YScale','linear');
end
set(gca,'FontName','Times New Roman')

% Axis labels and title
title('Entire area');
xlabel(xname, 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Count', 'FontWeight', 'bold', 'FontSize', 16);


if strcmp(choice, 'BdB_l')
xlabel(xname, 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Count', 'FontWeight', 'bold', 'FontSize', 18);
end


% Show mean and median
meanVal = mean(entireValsPos);
medVal  = median(entireValsPos);

hMed = xline(medVal, 'b--', 'LineWidth', 1.5);
hMed.DisplayName = sprintf('Median = %.2f', medVal);

hMean = xline(meanVal, 'k--', 'LineWidth', 1.5);
hMean.DisplayName = sprintf('Mean = %.2f', meanVal);




% Customize axes box
ax = gca;
ax.Box = 'on';  
ax.LineWidth = 0.5;  
ax.XColor = 'k';  
ax.YColor = 'k';  
ax.TickDir = 'in';  



% Add legend
legend([hMed, hMean], {hMed.DisplayName,hMean.DisplayName}, 'Location', 'northwest');

if strcmp(choice, 'BdB_l')
    xlim([2*10^-1, 2*10^8]);   
    ylim([0, 4.5*10^5]); 
    set(gca, 'YScale','linear');
    set(gca, 'FontSize',16);
    set(gca,'TickLength',[0.02, 0.002])
    ax.LineWidth = 1.5;  
    legend('hide');
    title('');
    xticks([1 10 100 1000 10000 100000 1000000 10000000 100000000])
    yticks([0 1*10^5 2*10^5 3*10^5 4*10^5])
end

hold off;

%  (b) Closed area
nexttile;% Second histogram
hold on;

% Compute histogram data
[counts_closed, binEdges] = histcounts(closedRegionValsPos, edges);
cumulative_closed = cumsum(counts_closed) / sum(counts_closed);

for i = 1:numBins
    b = bar(binCenters(i), counts_closed(i), 'FaceColor', colormapData(colorIndices(i), :), ...
        'EdgeColor', 'k', 'BarWidth', binWidths(i));
    if strcmp(choice, 'BdB_l')
        b.LineWidth = 1;
    end
end

% Set X-axis scale
if strcmp(choice, 'B') %|| strcmp(choice, 'BdB_l')
    set(gca, 'XScale','linear');
else
    set(gca, 'XScale','log');
    xlim([2*10^-7, 10^10]);
end

% Plot cumulative line
yyaxis right;
set(gca, 'Ycolor','b');
plot(binCenters, cumulative_closed, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Cumulative');
ylabel('Cumulative Probability','color','b', 'FontWeight', 'bold');
yyaxis left;

% Set Y-axis scale
set(gca, 'YScale','log');
    ylim([1, 10^6]);
if strcmp(choice, 'BdB_l')
    set(gca, 'YScale','linear');
end
set(gca,'FontName','Times New Roman')

title('Closed area');
xlabel(xname, 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Count', 'FontWeight', 'bold', 'FontSize', 16);


if strcmp(choice, 'BdB_l')
xlabel(xname, 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Count', 'FontWeight', 'bold', 'FontSize', 18);
end


% Display mean, median, and top 1%
meanVal = mean(closedRegionValsPos);
medVal  = median(closedRegionValsPos);
top1Val = prctile(closedRegionValsPos, 99); % Top 1%        

percentile_mean = 100 * sum(closedRegionValsPos < meanVal) / numel(closedRegionValsPos);
topPercent_mean = 100 - percentile_mean;

hMed = xline(medVal, 'b--', 'LineWidth', 1.5);
hMed.DisplayName = sprintf('Median = %.2f', medVal);

hMean = xline(meanVal, 'k--', 'LineWidth', 1.5);
hMean.DisplayName = sprintf('Mean = %.2f (Top %.1f%%)', meanVal, topPercent_mean);

hTop1 = xline(top1Val, 'r--', 'LineWidth', 1.5);  % 빨간 점선으로 표시
hTop1.DisplayName = sprintf('Top 1%% = %.2f', top1Val);


% Customize axes box
ax = gca;
ax.Box = 'on';  
ax.LineWidth = 0.5; 
ax.XColor = 'k';
ax.YColor = 'k'; 
ax.TickDir = 'in';  

% Add legend
legend([hMed, hMean, hTop1], ...
    {hMed.DisplayName, hMean.DisplayName, hTop1.DisplayName}, ...
    'Location', 'northwest');

if strcmp(choice, 'BdB_l')
    xlim([2*10^-1, 2*10^8]);   
    ylim([0, 10^4]); 
    set(gca, 'YScale','linear');
    set(gca, 'FontSize',16);
    set(gca,'TickLength',[0.02, 0.002])
    ax.LineWidth = 1.5;  % 테두리 두께 설정
    legend('hide');
    title('');
    xticks([1 10 100 1000 10000 100000 1000000 10000000 100000000])
    yticks([0 2*10^3 4*10^3 6*10^3 8*10^3 10*10^3])
end


% Add colorbar (log scale between 1 and 10^7)
colormap(turbo);
caxis([1, 10^7]);  
c = colorbar;
c.Ticks = logspace(0, 7, 8);
c.Label.String = xname;
c.Label.Interpreter = 'latex'; 
c.LineWidth = 1.5;

set(gca, 'ColorScale', 'log');  


hold off;



sgtitle(sprintf('Histograms with Mean & Median (D=%g μm) ', cd));

%% histogram save
%
if strcmp(choice, 'BdB_l')

else
csvFile = sprintf(['histogram_closed_d_%g_rot_%g_', choice, '.csv'], cd,rotationAngle);
histData = table(binStart', binEnd', binCenters', counts_entire', cumulative_entire', counts_closed', cumulative_closed',...
    'VariableNames', {'Bin_Start', 'Bin_End', 'Bin_Center', 'Entire_Count','Entire_cumulative', 'Closed_Count','Closed_cumulative'});

writetable(histData, csvFile);
end
%}

%% (6) Visualization with imagesc (NaN values shown as white)
%
if strcmp(choice, 'B') || strcmp(choice, 'BdB')


Z_entire2D = Z2D;   % Use original Z2D

figure2 = figure;

set(figure2, 'Units', 'normalized', 'OuterPosition', [0,0,0.4,0.6]);
h4 = imagesc(xUni, yUni, Z_entire2D);
set(gca, 'YDir','normal');  
axis equal;      % Keep aspect ratio
pbaspect([1 1 1]);  % Fix X:Y to 1:1
caxis(color_limit);  % Use predefined color limit
colormap(turbo);
c=colorbar;
c.Label.String = xname;
c.LineWidth=1.2;
c.Label.Interpreter = 'latex'; 


set(gca, 'ColorScale', ll);  % linear or log scale
%title(sprintf('Inside of the beads (including boundary, D=%g μm)', cd));
fontsize(18,"points");
xlabel('x [mm]','fontweight','bold','fontsize',20);
ylabel('y [mm]','fontweight','bold','fontsize',20);
set(gca, 'FontWeight', 'bold');
set(gca,'FontName','Times New Roman')

% Treat NaN as transparent, set background to white
set(h4, 'AlphaData', ~isnan(Z_entire2D));
set(gca, 'Color', 'w');



hold on;

if strcmp(choice,'B')
    viscircles(circleCenters_rotated, r,'Color', "#D4D4D4", 'LineWidth', 1);
end

% Generate magnetic flux lines
startx = linspace(min(X(:)),max(X(:)),20);  % Seed points (x)
starty = repmat(min(Y(:)), size(startx));  % Seed points (y)
verts = stream2(xUni,yUni, Bx2D, By2D, startx, starty, [0.1,20000]); 
hstream=streamline(verts);
set(hstream,'color','k','LineWidth',0.15);

ax = gca;
ax.Box = 'on';  
ax.LineWidth = 0.5;  
ax.XColor = 'k';  
ax.YColor = 'k';  
ax.TickDir = 'in';  
set(gca,'linewidth',1.2)

hold off;

%%
Z_closeed2D = nan(size(Z2D));
Z_closeed2D(closedRegionMask) = Z2D(closedRegionMask);


figure3 = figure;
hold on;
%viscircles(circleCenters_rotated, r,'Color', "#D4D4D4", 'LineWidth', 1.5);
set(figure3, 'Units', 'normalized', 'OuterPosition', [0,0,0.4,0.6]);

h5 = imagesc(xUni, yUni, Z_closeed2D);
set(gca, 'YDir','normal');  
axis equal;
pbaspect([1 1 1]); 
caxis(color_limit);  
colormap(turbo);
c=colorbar;
c.Label.String = xname;
c.LineWidth=1.2;
c.Label.Interpreter = 'latex'; 

set(gca, 'ColorScale', ll); % linear or log



%title(sprintf('Inside of the beads (including boundary, D=%g μm)', cd));
fontsize(18,"points");
xlabel('x [mm]','fontweight','bold','fontsize',20);
ylabel('y [mm]','fontweight','bold','fontsize',20);
set(gca, 'FontWeight', 'bold');
set(gca,'FontName','Times New Roman')

% Treat NaN as transparent, set background to white
set(h5, 'AlphaData', ~isnan(Z_closeed2D));
set(gca, 'Color', 'w');

ax = gca;
ax.Box = 'on';  
ax.LineWidth = 0.5;  
ax.XColor = 'k';  
ax.YColor = 'k';  
ax.TickDir = 'in';  
set(gca,'linewidth',1.2)
hold off;

end
%}
%% save
%

histFile = sprintf(['histogram  d_%g rot_%g ' ,choice,'.jpg'], cd,rotationAngle);
print(figure1, '-dpng', ['-r' num2str(300)], histFile);

%
if strcmp(choice, 'B') || strcmp(choice, 'BdB')
inFile = sprintf(['map in d_%g rot_%g ' ,choice,'.jpg'], cd,rotationAngle);
print(figure2, '-dpng', ['-r' num2str(300)], inFile);
inFile = sprintf(['map all d_%g rot_%g ' ,choice,'.jpg'], cd,rotationAngle);
print(figure3, '-dpng', ['-r' num2str(300)], inFile);
end
%}
%}