function p=testTaskRelated(falign,group_id)

%falign: nimage x ntrial
a=unique(group_id);
ss=size(falign);
fbin=zeros(length(a),ss(2));
group=zeros(size(fbin));
for i=1:length(a)
    id=(group_id==a(i));
    fbin(i,:)=mean(falign(id,:));
    group(i,:)=i;
end
    
p = kruskalwallis(fbin,group)