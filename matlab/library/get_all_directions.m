% Find all directions given a max distance
% This code can be optimized.
function [directions,normalized_directions] = get_all_directions(radius, dim)

if (dim == 3)
    % Brute force
    directions = [];
    for dx = -floor(radius):ceil(radius)
        for dy = -floor(radius):ceil(radius)
            for dz = -floor(radius):ceil(radius)
                
                if (dx == 0) && (dy == 0) && (dz == 0)
                    continue;
                end
                
                if sqrt( dx*dx +dy*dy + dz*dz ) <= radius
                    dir = [dx dy dz];
                    directions = keep_smallest(directions, dir, dim);
                end
            end
        end
    end
    
    normalized_directions = directions ./ repmat(sqrt(directions(:,1).^2 + directions(:,2).^2 + directions(:,3).^2), [1 3]);
    
elseif (dim == 2)
    
    directions = [];
    for dx = -floor(radius):ceil(radius)
        for dy = -floor(radius):ceil(radius)
            
            if (dx == 0) && (dy == 0)
                continue;
            end
            
            if sqrt( dx*dx +dy*dy) <= radius
                dir = [dx dy];
                directions = keep_smallest(directions, dir, dim);
            end
        end
    end
    
    normalized_directions = directions ./ repmat(sqrt(directions(:,1).^2 + directions(:,2).^2), [1 2]);
else
    error('Only 2 and 3 dimensions supported \n');
end

directions = int32(directions);

function directions = keep_smallest(directions, new, dim)


duplicate = 0;
for i = 1:size(directions,1)
    ndi = directions(i,:)./repmat(norm(directions(i,:)),[1 dim]);
    ndj = new ./repmat(norm(new),[1 dim]);
    
    % Same direction?
    if norm(ndi - ndj) < 1e-8
        duplicate = i;
        break;
    end
end

if (duplicate ~= 0)
    if (norm(new) < norm(directions(duplicate,:)))
        directions(duplicate,:) = new;
    end
else
    directions = [directions; new];
end