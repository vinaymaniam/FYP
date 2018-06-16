function [outputArg1,outputArg2] = plotme(x,y,c,type)
    lc = sprintf('%sx-',c);
    if type == 1
        if strcmp(c,'mix')
            plot(x,y,'x-','MarkerSize',8,'LineWidth',2)
        else
            plot(x,y,lc,'MarkerSize',8,'MarkerEdgeColor',c,'MarkerFaceColor',c,...
                'LineWidth',2)
        end
    elseif type == 2
        if strcmp(c,'mix')
            semilogx(x,y,'x-','MarkerSize',8,'LineWidth',2)
        else
            semilogx(x,y,lc,'MarkerSize',8,'MarkerEdgeColor',c,'MarkerFaceColor',c,...
                'LineWidth',2)
        end
end

