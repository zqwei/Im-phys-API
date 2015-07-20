function Iout=imshift(Iin,shift1,shift2)
ss=size(Iin);
Iout=zeros(ss);
idx1=(1:ss(1))-shift1;
idx2=(1:ss(2))-shift2;
idx1(idx1>ss(1))=[];
idx1(idx1<1)=[];
idx2(idx2>ss(2))=[];
idx2(idx2<1)=[];
Iout(idx1,idx2)=Iin(idx1+shift1,idx2+shift2);
    