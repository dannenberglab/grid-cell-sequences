function [grid_maps] = get_gridMaps(nmbr_neurons)

% This function computes all possible arrangements of cells in a "grid map"
% that are consistent with a trajectory code in 2D space by cell sequences
% given the number of neurons (nmbr_neurons) that are available to create
% the grid map in a rectangular environment. Input is the number of neurons
% that participate in the grid map. Output is a cell that contains all
% possible grid maps. These possibilities include reflections or 90 degree rotations and
% reflections of each other.

grid_maps = cell(1,1); % initiate cell with unknown number of elements

divisors = findDivisors(nmbr_neurons);

% first, try disjoint repeats across rows
for k = 1:numel(divisors)
    k2 = nmbr_neurons/divisors(k); % second factor so that divisors(k) * k2 = nmbr_neurons; k2 = number of rows that can be filled with disjoint neuron repeats
    for rows = 1:k2
        sequence{rows} = (rows-1)*divisors(k)+1:divisors(k)*rows;
    end
    if rows == 1
        % try translation (from 4 to nmbr_neurons-1); NOTE: the numbers
        % here indicate the INDEX of the cell in the translated repeat that
        % lies to the upper right of the first index of the repeat in the original row.
        for t = 4:nmbr_neurons-1
            if (t-1)*2 == nmbr_neurons || (t-1)*2 == nmbr_neurons + 1 || (t-1)*2 == nmbr_neurons + 2 % in these cases, the third row is translated in a way so that it creates an ambiguous code
                translation{t-3} = [];
                continue
            end
            translation{t-3} = [sequence{1}(t:end) sequence{1}(1:t-1)];
            % now continue in the same way but introduce a shift every
            % second row to plot the grid map aligned to a vertical wall
            temp = translation{t-3};
            for grid_rows = 2:20
                if iseven(grid_rows)
                    more_translations{grid_rows} = [temp(t-1:end) temp(1:t-2)];
                else
                    more_translations{grid_rows} = [temp(t:end) temp(1:t-1)];
                end
                temp = more_translations{grid_rows};
            end
            grid_maps{1,t-3} = repmat([sequence{1}; translation{t-3}; more_translations{2};more_translations{3};more_translations{4};more_translations{5};more_translations{6};more_translations{7};more_translations{8};more_translations{9};more_translations{10};more_translations{11};more_translations{12};more_translations{13};more_translations{14}; more_translations{15}; more_translations{16}; more_translations{17}; more_translations{18}; more_translations{19}; more_translations{20}],1,1);
            % plot grid maps
            plot_grid(grid_maps{1,t-3})
        end
    end
    
    if rows == 2
        % try translation (from 3 to nmbr_neurons-1)
        for t = 3:numel(sequence{1})-1
            translation{t-2} = [sequence{1}(t:end) sequence{1}(1:t-1)];
            % now continue
            temp = sequence{2};
            more_translations{3} = translation{t-2}; %to occupy cell for next for loop
            for grid_rows = 4:20
                more_translations{grid_rows} = [temp(t:end) temp(1:t-1)];
                temp = more_translations{grid_rows-1};
            end
            grid_maps{1,t-2} = repmat([sequence{1}; sequence{2}; translation{t-2};more_translations{4};more_translations{5};more_translations{6};more_translations{7};more_translations{8};more_translations{9};more_translations{10};more_translations{11};more_translations{12};more_translations{13};more_translations{14}; more_translations{15}; more_translations{16}; more_translations{17}; more_translations{18}; more_translations{19}; more_translations{20}],1,1);
            % plot grid maps
            plot_grid(grid_maps{1,t-2})
        end
    end
    
    if rows == 3
        % all translations of the first row work for the fourth row
        for t = 1:numel(sequence{1})
            translation{t} = [sequence{1}(t:end) sequence{1}(1:t-1)];
            % now continue in the same way but introduce a shift every
            % second row to plot the grid map aligned to a vertical wall
            temp = sequence{2};
            more_translations{3} = sequence{3}; %to occupy cell for next for loop
            more_translations{4} = translation{t}; %to occupy cell for next for loop
            for grid_rows = 5:20
                if ~iseven(grid_rows)
                    if t == 1
                        more_translations{grid_rows} = [temp(end) temp(1:end-1)];
                    elseif t > 1
                        more_translations{grid_rows} = [temp(t-1:end) temp(1:t-2)];
                    end
                elseif iseven(grid_rows)
                    more_translations{grid_rows} = [temp(t:end) temp(1:t-1)];
                end
                temp = more_translations{grid_rows-2};
            end
            grid_maps{1,t} = repmat([sequence{1}; sequence{2}; sequence{3}; translation{t};more_translations{5};more_translations{6};more_translations{7};more_translations{8};more_translations{9};more_translations{10};more_translations{11};more_translations{12};more_translations{13};more_translations{14}; more_translations{15}; more_translations{16}; more_translations{17}; more_translations{18}; more_translations{19}; more_translations{20}],1,1);
            % plot grid maps
            plot_grid(grid_maps{1,t})
        end
    end
    
    if rows == 4
        % all translations of the first row work for the fifth row
        for t = 1:numel(sequence{1})
            translation{t} = [sequence{1}(t:end) sequence{1}(1:t-1)]; % first translation in row 5 (translation of first row)
            % now continue
            temp = sequence{2};
            more_translations{3} = sequence{3}; %to occupy cell for next for loop
            more_translations{4} = sequence{4};
            more_translations{5} = translation{t};
            for grid_rows = 6:20
                more_translations{grid_rows} = [temp(t:end) temp(1:t-1)];
                temp = more_translations{grid_rows-3};
            end
            grid_maps{1,t} = repmat([sequence{1}; sequence{2}; sequence{3}; sequence{4}; translation{t}; more_translations{6};more_translations{7};more_translations{8};more_translations{9};more_translations{10};more_translations{11};more_translations{12};more_translations{13};more_translations{14}; more_translations{15}; more_translations{16}; more_translations{17}; more_translations{18}; more_translations{19}; more_translations{20}],1,1);
            % plot grid maps
            plot_grid(grid_maps{1,t})
        end
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper functions
    function divisors = findDivisors(number)
        divisors = [];
        for i = 3:number % smallest possible repeat length is 3
            if mod(number, i) == 0
                divisors = [divisors, i];
            end
        end
    end

    function [] = plot_grid(sequences)
        % Define parameters
        radius = 0.6; % radius of the circles
        spacing = 1.2; % spacing between centers of adjacent circles
        numRows = 20; % number of rows in the lattice
        numCols = 20; % number of columns in the lattice
        
        % Compute positions of circle centers
        x = [];
        y = [];
        for i = 1:numRows
            for j = 1:numCols
                if mod(i,2) == 1 % odd row
                    x(end+1) = (j-1)*spacing;
                    y(end+1) = (i-1)*spacing*sqrt(3)/2;
                else % even row
                    x(end+1) = (j-1)*spacing + spacing/2;
                    y(end+1) = (i-1)*spacing*sqrt(3)/2;
                end
            end
        end
        
        for row = 1:numRows
            for column = 1:numCols
                c = (row -1)*numCols + column;
                %                 color{c} = colors{sequences(row,mod(column,size(sequences,2))+1)};
                gridcell{c} = num2str(sequences(row,mod(column,size(sequences,2))+1));
            end
        end
        
        % Plot circles
        figure
        for i = 1:length(x)
            rectangle('Position', [x(i)-radius, y(i)-radius, 2*radius, 2*radius], 'Curvature', [1 1])
            % color one cell for better visualization
            if str2num(gridcell{i}) == 1
                rectangle('Position', [x(i)-radius, y(i)-radius, 2*radius, 2*radius], 'Curvature', [1 1], 'FaceColor', [0.5 0.5 0.5])
            end
            % Plot grid cell numbers within the grid
            text(x(i), y(i), gridcell{i}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12);
            % Set axis limits and aspect ratio
            axis([min(x)-radius max(x)+radius min(y)-radius max(y)+radius])
            axis equal
        end
        
    end

end
