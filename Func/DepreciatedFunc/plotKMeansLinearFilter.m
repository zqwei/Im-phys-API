function plotKMeansLinearFilter(dataStruct, k, params)
    neuronNum = size(dataStruct,1);
    timePointNum = size(dataStruct(1).unit_yes_trial_linear,2);
    shapeMat = zeros(neuronNum,timePointNum*2);
    for ii=1:neuronNum
      shapeMat(ii,1:timePointNum) = mean(dataStruct(ii).unit_yes_trial_linear);
      shapeMat(ii,(timePointNum+1):end) = mean(dataStruct(ii).unit_no_trial_linear);
    end
    shapeMat(:,end) = [];
    shapeMat(:,timePointNum) = [];
    timePointNum = timePointNum-1;

    % shapeMat = shapeMat./repmat(max(shapeMat,[],2),1,timePointNum*2);

    pcNum = k-1;
    [pc,score,~,~] = pca(shapeMat,'Algorithm','eig','Economy',false);

    clusterNum = k;
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
      plot(params.timeSeries(1:timePointNum), centroids(1:timePointNum,ii), '-', 'linewid', 1.0, 'color', [0.7 0 0]);
      hold on, plot(params.timeSeries(1:timePointNum), centroids((timePointNum+1):end,ii), '-', 'linewid', 1.0, 'color', [0 0 0.7])
      set(gca,'xlim',[params.timeSeries(1) params.timeSeries(end-1)])
      gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
      title(['% neurons ' num2str(100*neuronsPerCluster(ii)/neuronNum)])
      box off
    end

end