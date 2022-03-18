function G = complete_graph(nvert)
%COMPLETE_GRAPH Generates a Matlab graph object that represents a complete
%graph with a given numer ob vertices.
%   This function genenerates a Matlab graph object that contains a
%   ccomplete graph with the given amount of vertices.
%
%   Arguments:
%       nvert    -> Number of vertices

G = graph(ones(nvert) - eye(nvert));
end
