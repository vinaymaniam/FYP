function [out] = flipimg(img,mode)
    out = img;
    switch mode
        case 2
            out = fliplr(img);            
        case 3
            out = flip(img);
        case 4
            out = flip(img);
            out = fliplr(img);
    end
end

