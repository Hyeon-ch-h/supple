clear; clc;

%% (0) Read data


% Ask user for the input filename (without .csv)
filename = "(filename)";
choice = ['BdB'];     % 'B','BdB','dB2','Sm','BdB_l'  / choose one


%% (1) Seprate Data
% Define full file paths
input_file = filename + ".csv";            % Input CSV file
% Read the CSV file (skipping header)
loaded= readtable(input_file);

% Z depending on choice
if strcmp(choice, 'B')
    Z = loaded{:, 3}; 
    xname = 'B [T]';  % x label on histogram
    color_limit=[0, 0.7];
    ll='linear';
elseif strcmp(choice, 'BdB') || strcmp(choice, 'BdB_l')
    Z = loaded{:, 3}/2; 
    xname = '(B·∇)B [T^2/m]';
    color_limit=[1, 10^7];
    ll='log';
elseif strcmp(choice, 'dB2')
    Z = loaded{:, 3}; 
    xname = '∇B^2 [T^2/m]';
    color_limit=[1, 10^7];
    ll='log';
elseif strcmp(choice, 'Sm')
    Z = loaded{:, 3}/2/(4*pi*10^(-7)); 
    xname = 'S_m [TA/m^2]';
    color_limit=[1, 10^7];
    ll='log';
else
    error('Invalid choice. Please select a valid option.');
end



%%
% --- Set histogram bins ---
numBins = 101;

if strcmp(choice, 'B') %|| strcmp(choice, 'BdB_l')
    minV=0;
    maxV=2;
    edges = linspace(minV, maxV, numBins);  % Linear binning
    % Arithmetic mean of linear bin edges
    binCenters = (edges(1:end-1) + edges(2:end)) / 2;
else
    minV=10^-6;
    maxV=10^11;
    edges = logspace(log10(minV), log10(maxV), numBins);  % Logarithmic binning
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
binCenters_fixed = max(binCenters, 1);   % Replace values < 1 with 1 for log transform
log_binCenters = log10(binCenters_fixed);  

log_min = log10(1);      
log_max = log10(10^7);   
log_binCenters = (log_binCenters - log_min) / (log_max - log_min);  

% Adjust color indices to prevent exceeding colormap range
colorIndices = round(1 + log_binCenters * (numBins - 1));  
colorIndices = min(colorIndices, numBins);  



%% Plot
figure1 = figure('Name','Histograms: Entire / Excluded / Remain','NumberTitle','off');

set(gcf, 'Units', 'normalized', 'OuterPosition', [0,0,0.4,0.8]);

if strcmp(choice, 'BdB_l')
set(gcf, 'Units', 'normalized', 'OuterPosition', [0,0,0.35,0.5]);
end

hold on;


[counts_entire, binEdges] = histcounts(Z, edges);

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

if strcmp(choice, 'B') %|| strcmp(choice, 'BdB_l')
    set(gca, 'XScale','linear');
else
    set(gca, 'XScale','log');
    xlim([2*10^-7, 10^10]);
end

yyaxis right;
set(gca, 'YColor', 'b');
plot(binCenters, cumulative_entire, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Cumulative');
ylabel('Cumulative Probability','color','b');


yyaxis left;
set(gca, 'YScale','log');
    ylim([1, 10^6]);
if strcmp(choice, 'BdB_l')
    set(gca, 'YScale','linear');
end

title('Closed area');
xlabel(xname, 'FontWeight', 'bold', 'FontSize', 14);
ylabel('Count', 'FontWeight', 'bold', 'FontSize', 14);


if strcmp(choice, 'BdB_l')
xlabel(xname, 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Count', 'FontWeight', 'bold', 'FontSize', 16);
end


meanVal = mean(Z);
medVal  = median(Z);
top1Val = prctile(Z, 99); 

percentile_mean = 100 * sum(Z < meanVal) / numel(Z);
topPercent_mean = 100 - percentile_mean;


hMed = xline(medVal, 'b--', 'LineWidth', 1.5);
hMed.DisplayName = sprintf('Median = %.2f', medVal);

hMean = xline(meanVal, 'k--', 'LineWidth', 1.5);
hMean.DisplayName = sprintf('Mean = %.2f (Top %.1f%%)', meanVal, topPercent_mean);

hTop1 = xline(top1Val, 'r--', 'LineWidth', 1.5);  
hTop1.DisplayName = sprintf('Top 1%% = %.2f', top1Val);


ax = gca;
ax.Box = 'on';  
ax.LineWidth = 0.5; 
ax.XColor = 'k';  
ax.YColor = 'k';  
ax.TickDir = 'in'; 


legend([hMed, hMean, hTop1], ...
    {hMed.DisplayName, hMean.DisplayName, hTop1.DisplayName}, ...
    'Location', 'northwest');

if strcmp(choice, 'BdB_l')
    xlim([2*10^-1, 2*10^8]);   
    ylim([0, 10^4]); 
    set(gca, 'YScale','linear');
    set(gca, 'FontSize',14);
    set(gca,'TickLength',[0.02, 0.002])
    ax.LineWidth = 1.5;  % 테두리 두께 설정
    legend('hide');
    title('');
    xticks([1 10 100 1000 10000 100000 1000000 10000000 100000000])
    yticks([0 2*10^3 4*10^3 6*10^3 8*10^3 10*10^3])
end



colormap(turbo);
caxis([1, 10^7]);  
c = colorbar;
c.Label.String = xname;
set(gca, 'ColorScale', 'log');  


hold off;


%% histogram save
%
if strcmp(choice, 'BdB_l')

else
csvFile = sprintf(['%s histogram.csv'],filename);
histData = table(binStart', binEnd', binCenters', counts_entire', cumulative_entire', ...
    'VariableNames', {'Bin_Start', 'Bin_End', 'Bin_Center', 'Entire_Count', 'Cumulative'});

writetable(histData, csvFile);
end
%}

%% save
%

histFile = sprintf(['%s histogram ' ,choice,'.jpg'], filename);
print(figure1, '-dpng', ['-r' num2str(300)], histFile);

%}