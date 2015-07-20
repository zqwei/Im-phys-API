function ROI_list_out=make_neuropil_mask(ROI_list)    
    ss=[512,512];
    ss=ss(1:2);
    
    ROI_map=zeros(ss(1:2));
    for j=1:length(ROI_list)
        ROI_map(ROI_list(j).pixel_list)=1;
    end
%     figure;imagesc(ROI_map);
    
    ROI_map=imresize(ROI_map,ss(1:2));
    se = strel('disk', 5,0);
    ROI_map_dilate=imdilate(ROI_map,se);
    ROI_map_dilate=imfill(ROI_map_dilate,'holes');
    
    ROI_list_out=[];
    radius=40;
    for j=1:length(ROI_list)
        Neuropil_map=zeros(ss(1:2));      
        center=round(ROI_list(j).centerPos);
        
        for ii=-radius:radius
            for jj=-radius:radius
                if (sqrt(ii^2+jj^2)<radius)  && (center(1)+ii)>0  && (center(1)+ii)<ss(1)  && (center(2)+jj)>0   && (center(2)+jj)<(ss(2))
                Neuropil_map(center(1)+ii,center(2)+jj)=1;
                end
            end
        end
         
        Neuropil_map=Neuropil_map & (~ROI_map_dilate);
        Neuropil_map=imresize(Neuropil_map,ss);
%         figure;imagesc(Neuropil_map);
        temp=ROI_list(j);        
        temp.neuropil_list=find(Neuropil_map);
        
%         sub=im(temp.Neuropil_list,:);
        temp.fmean_neuropil=[];
        ROI_list_out=[ROI_list_out,temp];   
    end
    
    
end
%     im=readTifStack([name(1:end-4),'.tif']);
    