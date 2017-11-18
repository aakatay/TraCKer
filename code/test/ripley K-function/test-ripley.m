v = 1:10; N = 10; DX = repmat(v',1,N)-repmat(v,N,1);

%%
training = [mvnrnd([ 1  1],   eye(2), 100); ...
            mvnrnd([-1 -1], 2*eye(2), 100)];
group = [repmat(1,100,1); repmat(2,100,1)];
gscatter(training(:,1),training(:,2),group,'rb','+x');
legend('Training group 1', 'Training group 2');

%%
       % training data: two normal components
       training = [mvnrnd([ 1  1],   eye(2), 100); ...
                   mvnrnd([-1 -1], 2*eye(2), 100)];
       group = [ones(100,1); 2*ones(100,1)];
       gscatter(training(:,1),training(:,2),group);hold on;
 
       % some random sample data
       sample = unifrnd(-5, 5, 100, 2);
       % classify the sample using the nearest neighbor classification
       c = knnclassify(sample, training, group);
 
       gscatter(sample(:,1),sample(:,2),c,'mc'); hold on;
       c3 = knnclassify(sample, training, group, 3);
       gscatter(sample(:,1),sample(:,2),c3,'mc','o');

%%
       a =.9;
       mu = [1 -1]; Sigma = [.9 .4; .4 .3];
       mu = [1 -1]; Sigma = [.9 a; a .3];
       r = mvnrnd(mu, Sigma, 500);
       plot(r(:,1),r(:,2),'.');
       set(gca,'xlim',[-3 5],'ylim',[-3 1])
       
       
       %% kmeans
       
       
 
        X = [randn(20,2)+1.5*ones(20,2);randn(20,2)+ones(20,2); randn(20,2)-0.5*ones(20,2)];
        opts = statset('Display','final');
        [cidx, ctrs] = kmeans(X, 5, 'Distance','city', ...
                              'Replicates',5, 'Options',opts);
        plot(X(cidx==1,1),X(cidx==1,2),'r.', ...
             X(cidx==2,1),X(cidx==2,2),'y.', ...
             X(cidx==3,1),X(cidx==3,2),'b.', ctrs(:,1),ctrs(:,2),'kx');