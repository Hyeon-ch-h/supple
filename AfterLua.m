clear; clc;

% Define the base filename (without extension)
filename = "bead 250 45 magnet off";  % Change this to your filename if needed

% Define full file paths
input_file = filename + ".csv";             % Input CSV file
output_file = filename + "_grad.csv";         % Output CSV file with "_grad" suffix

% Read the CSV file (assuming a header exists)
data = readmatrix(input_file);

% Extract columns (assumed order: X, Y, Bx, By, ... )
X = data(:,1);  % X coordinates (in mm)
Y = data(:,2);  % Y coordinates (in mm)
Bx = data(:,3); % Bx component
By = data(:,4); % By component

% Compute magnetic field magnitude
Bmag = sqrt(Bx.^2 + By.^2);

% Compute Bmag^2
Bmag_squared = Bmag .^ 2;

% Get unique X and Y values (assumes a regular grid)
unique_x = unique(X);
unique_y = unique(Y);

% Automatically compute grid spacing dx and dy (in mm)
dx = (unique_x(end) - unique_x(1)) / (numel(unique_x)-1);
dy = (unique_y(end) - unique_y(1)) / (numel(unique_y)-1);

% Reshape data into 2D grids (assuming they are structured grid data)
nx = numel(unique_x);
ny = numel(unique_y);

% Reshape the vectors into matrices
Bmag2_grid = reshape(Bmag_squared, ny, nx);

% For gradient calculation, convert dx and dy to meters (since 1 m = 1000 mm)
dx_m = dx / 1000;
dy_m = dy / 1000;

% Compute the gradient of Bmag^2 on the grid using the converted spacing
[dB2_dx, dB2_dy] = gradient(Bmag2_grid, dx_m, dy_m);

% Compute gradient magnitude of B^2 (∇B^2) in T^2/m
GradB2_magnitude = sqrt(dB2_dx.^2 + dB2_dy.^2);

% Flatten data back into column vectors for output
GradB2_final = GradB2_magnitude(:);

% Create final dataset: columns are X (mm), Y (mm), B(T), ∇B^2 (T^2/m), Bx (T), By (T)
final_data = [X, Y, Bx, By, Bmag, GradB2_final];

% Write header and data to a new CSV file
header = ["X(mm)","Y(mm)", "Bx(T)", "By(T)","B(T)","∇B^2(T^2/m)"];
writematrix(header, output_file);  % Write header
writematrix(final_data, output_file, 'Delimiter', ',', 'WriteMode', 'append');

disp("B, ∇B^2, Bx, and By calculation complete! Data saved to: " + output_file);
