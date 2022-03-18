function G = triangle_graph(baselength)
%TRIANGLE_GRAPH Generates a Matlab graph object that represents a
%triangular graph.
%   This function generates a Matlab graph object that represents a
%   triangular graph. The generated graphs are always undirected and can be
%   configured only by their size, which is specified by giving the number
%   of vertices in the bottommost row.

%% Setup intermediate variables
% Calculate total number of agents
N = baselength * (baselength+1) / 2;

% Initial distributing algorithm for the first row
row_start = 1; % Index of the first node of this row
row_end   = 1; % Index of the last node of this row
row_no    = 1; % Number of the row as counted from the top

%% Main algorithm
A = zeros(N);

for i = 1:N
    % Check if node is on the bottom row, and if not, add two edges
    if row_no ~= baselength
        A(i, i+row_no)   = 1;
        A(i, i+row_no+1) = 1;
    end
    
    % Check if node is on the left edge of the triangle, if not, add two
    % edges
    if i ~= row_start
        A(i, i-row_no) = 1;
        A(i, i-1)      = 1;
    end
    
    % Check if node is on the right edge. If so, switch to the next row,
    % otherwise add two more edges
    if i ~= row_end
        A(i, i-row_no+1) = 1;
        A(i, i+1)        = 1;
    else
        row_no    = row_no  + 1;
        row_start = row_end + 1;
        row_end   = row_end + row_no;
    end
end

G = graph(A);
end
