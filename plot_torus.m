%This script plots the hexagonal arrangement of firing fields that results
%from a sequence code of trajectories using multiple fundatmental grid
%units consisting of seven cells each.
% The total number of cells is a multiple of 7*(4^n), where n is a non-negative integer (0,1,2,...).

% Requirements:
% function distinguishable_colors.m
% (https://github.com/cortex-lab/MATLAB-tools/blob/master/distinguishable_colors.m)
% MATLAB Image Processing Toolbox

n = 4; % a non-negative integer (0,1,2,...) that determins the number of rows and number of neurons;
nmbrRows = 2^n;%a positive integer (1,2,...); determines the circumference of the doughnut (how thick the doughnut is)
nmbrNeurons = 7*(4^n); % must be a multiple of 7*(4^n) with n = a non-negative integer (0,1,2,...)
nmbrColumns = 7*2^n; % same as neural sequence repeat length

% Define the parameters for the doughnut
R = (nmbrColumns+1)/2;        % Major radius
r = ((nmbrRows+1).*(sqrt(3)/2))/2;   % Minor radius (multiplied with cosd(30) = sqrt(3)/2)

theta = linspace(0, 2*pi-2*pi/nmbrRows, nmbrRows);    % Angular variable theta
phi = linspace(0, 2*pi-2*pi/nmbrColumns, nmbrColumns);      % Angular variable phi

% Generate the coordinates for the surface
[Theta, Phi_temp] = meshgrid(theta, phi); % phi is the big circle surrounding the hole in the center of the doughnut; theta is the small circle that determines how thick the doughnut is
% to create hexagonal lattice, shift every grid row (columns in Phi) by half the spacing between
% two cells
shift = (Phi_temp(2,1) - Phi_temp(1,1))/2;
Phi = Phi_temp;
for i = 2:size(Phi,2)
    Phi(:,i) = Phi(:,i-1) + shift;
end

X = (R + r*cos(Theta)) .* cos(Phi);
Y = (R + r*cos(Theta)) .* sin(Phi);
Z = r * sin(Theta);

% Create a figure
figure;
hold on;

neuralSequence = 1:7*2^n;
gridMap = nan(nmbrRows,numel(neuralSequence));
for i = 1:nmbrRows
    gridMap(i,:) = neuralSequence + neuralSequence(end)*(i-1);
end
gridMap = gridMap';

if nmbrNeurons <=28
    colorsAll = distinguishable_colors(28);
    stopRepeat = 28;
else
    colorsAll_temp = distinguishable_colors(numel(neuralSequence)*2,[1 1 1]);
    colorsAll = repmat(colorsAll_temp,nmbrRows,1);
    stopRepeat = 28*(log2(nmbrRows)-1)-nmbrRows;
end

colors = ones(numel(gridMap),3).*0.6;
if n > 0
    c = 0;
    for r = 1:stopRepeat
        for i = 1:numel(neuralSequence):numel(gridMap)
            c = c + 1;
            if c >= 14*nmbrRows
                break
            end
            colors(i+nmbrRows/2*(r),:) = colorsAll(c,:);
        end
    end
end

scatter3(X(:),Y(:),Z(:),125,colors,'filled')

% Set the view and axis labels
view(-60, 30);
xlabel('X');
ylabel('Y');
zlabel('Z');

% Set the title
title(sprintf('Neuronal Manifold (%d grid cells)', nmbrNeurons));

% Adjust the axis limits and aspect ratio
axis equal tight;

% Hide the grid lines
grid off;

% Remove the box around the plot
box off;

% Rotate the plot
rotate3d on;
