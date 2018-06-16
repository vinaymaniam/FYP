for stg = 1:2
    load(sprintf('Center2048_%istage.mat',stg));
    Center = single(Center);

    %% Kmeans with replicates to get better results
%     Center = Center';
    nctrds = round(sqrt(size(Center,1)));

    stream = RandStream('mlfg6331_64');
    options = statset('UseParallel',1,'UseSubstreams',1,...
    'Streams',stream);    
    [~, C] = kmeans(Center,nctrds,'Options',options,'MaxIter',300,...
                          'Display','final','Replicates',10);

    index = knnsearch(Center, C, 'K', 4*nctrds, 'Distance', 'euclidean');

    heirarchy = zeros(nctrds, 25, size(index,2)+1);
    heirarchy(:,:,1) = C;
    % ind = zeros(nctrds, size(index,2));
    for i = 1:size(index,2)
        for j = 1:nctrds
            heirarchy(j,:,i+1) = Center(index(j,i),:);
    %         ind(j,i) = index(j,i);
        end
    end

    % save('Heirarchy4096','heirarchy','index')
    save(sprintf('%iheir2048',stg),'heirarchy','index');
end