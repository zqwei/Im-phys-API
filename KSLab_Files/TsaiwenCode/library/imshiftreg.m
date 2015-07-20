function [Ireg,shift1,shift2,maxcorr,xcorrvalue]=imshiftreg(Imoving,Ifix,srange)
% quick and simple image registration by shifting integer pixel and
% maximize xcorr

ss=size(Ifix);

idx1=(srange+1):(ss(1)-srange);
idx2=(srange+1):(ss(2)-srange);
I1=Ifix(idx1,idx2);
normI1=I1(:)-mean(I1(:));
normI1=normI1/norm(normI1);

xcorrvalue=zeros(srange,srange);

range1=-srange:srange;
range2=-srange:srange;

for i=1:length(range1)   
    for j=1:length(range2)
        I2=Imoving(idx1+range1(i),idx2+range2(j));
        normI2=I2(:)-mean(I2(:));
        normI2=normI2/norm(normI2);
        xcorrvalue(i,j)=normI2'*normI1;
    end
end
% figure;imagesc(xcorrvalue);
[maxcorr,idx]=max(xcorrvalue(:));
[i,j]=ind2sub(size(xcorrvalue),idx);
shift1=i-srange-1;
shift2=j-srange-1;

%compensate
Ireg=zeros(ss);
idx1=(1:ss(1))-shift1;
idx2=(1:ss(2))-shift2;
idx1(idx1>ss(1))=[];
idx1(idx1<1)=[];
idx2(idx2>ss(2))=[];
idx2(idx2<1)=[];
Ireg(idx1,idx2)=Imoving(idx1+shift1,idx2+shift2);
