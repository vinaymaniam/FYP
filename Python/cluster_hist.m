load DX_and_DY/DX_all


for n = [256,1024,4096,16384]
    fprintf('n=%i\n',n)
    load(sprintf('data_files/pyCenter%i',n));

    centcounts = zeros(1,n);
    if n == 16384
        for i = 1:4096:16384
            ci = Center(:,i:i+4095);
            di = pdist2(X',ci');
            d = [d;di];
        end
    else
        d = pdist2(X',Center');
    end
    [~,idx] = min(d,[],2);
    for i = 1:length(idx)
        centcounts(idx(i)) = centcounts(idx(i)) + 1;
    end
    figure;
    ax = gca;
    ax.XLim = [0,500];
    histogram(ax,centcounts,100)
end
