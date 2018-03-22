for i = [256 1024 2048 4096]
    load(sprintf('pyMap%i',i));
    if iscell(Map)
        Map = cellfun(@(x)reshape(x,[],1),Map,'un',0);
        Map = cell2mat(Map);
        save(sprintf('pyMap%i',i),'Map');
    end
end