function G = grid_graph(width, height)
%GRID_GRAPH Summary of this function goes here
%   Detailed explanation goes here

if nargin <= 1
    height = width;
end

%% Generate graph structure
A = zeros(width*height);
for i = 1:width*height
    x = mod(i-1, width) + 1;
    y = floor((i-1)/width) + 1;
    
    if x > 1
        A(i,i-1) = 1;
    end
    if x < width
        A(i,i+1) = 1;
    end
    if y > 1
        A(i,i-width) = 1;
    end
    if y < height
        A(i,i+width) = 1;
    end
end

G = graph(A);
end
