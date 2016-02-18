function map=np_map(data)
    ss=size(data);
    data=double(data);
    neighbor_avg=data;
    H=[0.5,1,0.5;1,0,1;0.5,1,0.5];
    for i=1:ss(3)
        %neighbor_avg(:,:,i)=imfilter(neighbor_avg(:,:,i),H);
        neighbor_avg(:,:,i)=filter2(H,neighbor_avg(:,:,i));
    end
    data=reshape(data,[],ss(3));
    neighbor_avg=reshape(neighbor_avg,[],ss(3));
    
    neighbor_avg=sub_linear(neighbor_avg);
    data=sub_linear(data);
    
    map=zeros(size(data,1),1);
    for i=1:size(data,1)
        temp=neighbor_avg(i,:);
        temp=temp/norm(temp);
        map(i)=temp*data(i,:)'/norm(data(i,:));
    end
    map=reshape(map,ss(1:2));
end
    
    
    
