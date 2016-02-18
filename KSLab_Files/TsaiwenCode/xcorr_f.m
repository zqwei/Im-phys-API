function fmean_sub=xcorr_f(fmean)
n=floor(size(fmean,1)/3);
fmean=fmean(1:n*3,:);
ss=size(fmean);
fmean=reshape(fmean,3,[],ss(2));
fmean=mean(fmean);
fmean=reshape(fmean,[],ss(2));
b = hann(20);
b=b/sum(b);
fmean_sub=zeros(size(fmean));
for k=1:size(fmean,2)
    base=filtfilt(b,1,double(fmean(:,k)));
    fmean_sub(:,k)=fmean(:,k)-base;
    fmean_sub(:,k)=fmean_sub(:,k)-mean(fmean_sub(:,k));
end

d=pdist(fmean_sub','correlation');
e=squareform(d);
figure;imagesc(1-e)