function [y] = affine(x,k,flp,direction)
    if direction == 1
        k = 4-k;
        y = x;
        y = rot90(y,k);
        if flp > 0
            y = flip(y);
        end            
    else  
        if flp > 0
            y = flip(x);
        else
            y = x;
        end    
        y = rot90(y,k);
    end
end

