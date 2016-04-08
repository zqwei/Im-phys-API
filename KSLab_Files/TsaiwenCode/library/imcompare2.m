function imcompare2(I1,I2)
ss=size(I1);
hfig=figure('position',[100,100,1200,600],'WindowButtonMotionFcn',@MouseMotionFcn);
hax1=subplot(1,2,1);hold on;
if size(I1,3)==1
    him1=imagesc(I1);
else
    him1=image(I1);
end
hline1x=plot([0,ss(2)],[0,0],'-w');
hline1y=plot([0,0],[0,ss(1)],'-w');
axis image;set(gca,'ydir','reverse');

hax2=subplot(1,2,2);hold on;
if size(I2,3)==1
    him2=imagesc(I2);
else
    him2=image(I2);
end
hline2x=plot([0,ss(2)],[0,0],'-w');
hline2y=plot([0,0],[0,ss(1)],'-w');
axis image;set(gca,'ydir','reverse');
    function MouseMotionFcn(varargin)
        pt1=get(get(him1,'parent'),'currentPoint');
        pt2=get(get(him2,'parent'),'currentPoint');
        pt1=pt1(1,1:2);
        pt2=pt2(1,1:2);
        if pt1(1)<=ss(2) && pt1(1)>0 && pt1(2)<=ss(1) && pt1(2)>0
            pt=pt1;
        elseif pt2(1)<=ss(2) && pt2(1)>0 && pt2(2)<=ss(1) && pt2(2)>0
            pt=pt2;
        else
            pt=[];
        end     
        if ~isempty(pt);
            set(hline1x,'ydata',pt(2)*[1,1]);
            set(hline1y,'xdata',pt(1)*[1,1]);
            set(hline2x,'ydata',pt(2)*[1,1]);
            set(hline2y,'xdata',pt(1)*[1,1]);
        end
                    
    end
end