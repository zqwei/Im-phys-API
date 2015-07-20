function plot_multi_roi(t,dff1,dff2)
crange=[-1,1];
figure;
subplot(1,4,1);imagesc(t,[1,size(dff1,2)],dff1');set(gca,'clim',crange);hold on
plot([-2.6,-2.6],[0,size(dff1,2)],'k:','linewidth',2);
plot([-1.4,-1.4],[0,size(dff1,2)],'k:','linewidth',2);
plot([0,0],[0,size(dff1,2)],'k:','linewidth',2);
xlim([min(t),max(t)]);ylim([1,size(dff1,2)]);

subplot(1,4,2);imagesc(t,[1,size(dff2,2)],dff2');set(gca,'clim',crange);hold on
plot([-2.6,-2.6],[0,size(dff2,2)],'k:','linewidth',2);
plot([-1.4,-1.4],[0,size(dff2,2)],'k:','linewidth',2);
plot([0,0],[0,size(dff2,2)],'k:','linewidth',2);
xlim([min(t),max(t)]);ylim([1,size(dff2,2)]);

subplot(1,4,3);
h1=plot(t,squeeze(mean(dff1,2)),'color','b','linewidth',2);hold on;
h2=plot(t,squeeze(mean(dff2,2)),'color','r','linewidth',2);
plot([-2.6,-2.6],[-0.2,0.5],'k:','linewidth',2);
plot([-1.4,-1.4],[-0.2,0.5],'k:','linewidth',2);
plot([0,0],[-0.2,0.5],'k:','linewidth',2);box off;
xlim([min(t),max(t)]);
% legend([h1,h2],{'correct right','correct left'});

