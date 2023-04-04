% Solves maze with BFS
function [path,success,M,mid_finish] = solve_maze_ITER(img,start,mid_finish,M,v)
    % Solves a maze using Breadth First Search (BFS) algorithm
    
    % Input:
    % img: SAR image for conversion into binary maze
    % start: starting point in the maze [x, y]
    % mid_finish: kx2 matrix of k intermediate destination points in the maze
    % M: matrix of reachable set of pixels from start (optional)
    % v: verbosity level (0 for no output, 1 for output when path fails, 2 for output when path is successful)
    
    % Output:
    % path: matrix of points forming the shortest path from start to finish (if one exists)
    % success: boolean value indicating whether a path was found or not
    % M: matrix of reachable set of pixels from start
    % mid_finish: matrix of intermediate destination points in the maze
    
    success = 1;
    if isempty(M) % If M is not provided, create it
        maze = img > 0; % Convert image to binary

        %% Init BFS
        n = numel(maze);
        Q = zeros(n, 2);            % Initialize queue with size n
        M = zeros([size(maze) 2]);  % Initialize M with size [m x n x 2]
        front = 0;
        back = 1;

        [Q,M,front] = push(maze,M,Q,front,start,[0 0]);

        d = [0 1; 0 -1; 1 0; -1 0; 1 1; 1 -1;-1 1;-1 -1];   % Define the possible directions to move in the maze

        %% Run BFS
        while back <= front     % While there are still points to explore
            p = Q(back, :);     % Get the next point from the queue
            back = back + 1;
            for i = 1:4         % Check all 8 possible directions
                [Q,M,front] = push(maze,M,Q,front,p,d(i, :));   % Add point to queue and M if it is a valid path
            end
        end
    end

    %% Extracting path
    loop = 1;trigger=0;
    while loop<=size(mid_finish,1) && trigger == 0  % Check all intermediate destination points until a path is found   
        path = mid_finish(loop,:);
        while true
            q = path(end, :);
            if q(1)==0 || q(2)==0   % If the path fails, set success to 0 and try again with a new intermediate destination
                success=0;
                if v == 1           % Print a message if v is set to 1
                    [a,b]=find(M(:,:,1)>0);
                    [xa,Ia]=max(a);
                    ya=b(Ia);
                    mid_finish=[xa,ya];
                end
                break
            end
            p = reshape(M(q(1), q(2), :), 1, 2);    % Get the next point in the path from M
            path(end + 1, :) = p;                   % Add point to the path
            if isequal(p, start)                    % If the path reaches the starting point
                if v == 2
                    success=1;
                end
                trigger=1;
                mid_finish = mid_finish(loop,:);
                path = flip(path);
                break;
            end
        end
        loop=loop+1;
    end
end

function [Q,M,front]=push(maze,M,Q,front,p,d)
% PUSH updates the state of the search algorithm by adding a new point 
% to the queue and updating the state of the search. 
%
% Arguments:
%   - maze: A binary matrix that represents the maze. 
%   - M: Matrix of reachable set of pixels from start
%   - Q: A queue that stores the points to be visited by the algorithm.
%   - front: An integer for the position of the front of the queue.
%   - p: A vector that represents the current point in the search.
%   - d: A vector for the direction to be taken from the current point.
%
% Returns:
%   - Q: An updated queue that stores the points to be visited.
%   - M: Matrix of reachable set of pixels from start
%   - front: An integer for the position of the front of the queue.

    % Calculate the new point to be added to the queue
    q = p + d;

    % Check if the new point is a valid point to be visited by the algorithm
    if ~(q(1)==0 || q(2)==0) && ~(q(1)>= size(maze,1) || q(2)>= size(maze,2))
        if maze(q(1), q(2)) && M(q(1), q(2), 1) == 0
          % Update the front of the queue and add the new point to it
          front = front + 1;
          Q(front, :) = q;

          % Update the matrix M with the new path found
          M(q(1), q(2), :) = reshape(p, [1 1 2]);
        end
    end
end