function I_bp=backproj(SR,LR) 
    SR(1:2:end,1:2:end) = LR;
    I_bp = SR;
end

