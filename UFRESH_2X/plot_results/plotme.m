function [outputArg1,outputArg2] = plotme(x,y,c,type)
    lc = sprintf('%sx-',c);
    if type == 1
        plot(x,y,lc,'MarkerSize',8,'MarkerEdgeColor',c,'MarkerFaceColor',c,...
            'LineWidth',2)
    elseif type == 2
        semilogx(x,y,lc,'MarkerSize',8,'MarkerEdgeColor',c,'MarkerFaceColor',c,...
            'LineWidth',2)
end

