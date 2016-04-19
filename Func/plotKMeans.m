function plotKMeans(dataStruct)
    neuronNum = size(dataStruct,1);
    timePointNum = size(dataStruct(1).unit_yes_trial,2);
    shapeMat = zeros(neuronNum,timePointNum*2);
    for ii=1:neuronNum
      shapeMat(ii,1:timePointNum) = mean(dataStruct(ii).unit_yes_trial);
      shapeMat(ii,(timePointNum+1):end) = mean(dataStruct(ii).unit_no_trial);
    end
    shapeMat(:,end) = [];
    shapeMat(:,timePointNum) = [];
    timePointNum = timePointNum-1;

    % shapeMat = shapeMat./repmat(max(shapeMat,[],2),1,timePointNum*2);

    pcNum = 8;
    [pc,score,~,~] = pca(shapeMat);

    clusterNum = 9;
    [idx,C,~,~] = kmeans(score(:,1:pcNum),clusterNum);

    neuronsPerCluster = zeros(1,clusterNum);
    for ii=1:clusterNum
      neuronsPerCluster(ii) = sum(idx==ii);
    end

    centroids = pc(:,1:pcNum)*C';
    [~,ix] = sort(neuronsPerCluster,'descend');
    centroids = centroids(:,ix);
    neuronsPerCluster = neuronsPerCluster(ix);
    % Plot
    figure;
    squareNum = ceil(sqrt(clusterNum));
    for ii=1:clusterNum
      d = diff(centroids(:,ii)); d(timePointNum) = [];
      if mean(d)<0;
        centroids(:,ii) = -centroids(:,ii);
      end
      subplot(squareNum,squareNum,ii);
      plot(centroids(1:timePointNum,ii));
      hold on, plot(centroids((timePointNum+1):end,ii),'r')
      set(gca,'xlim',[0 timePointNum])
      title(['Perc. neurons ' num2str(100*neuronsPerCluster(ii)/neuronNum)])
    end

end