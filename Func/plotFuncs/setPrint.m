function setPrint(width, height, fname, printType)

%     h1 = suptitle(fname(9:end));
%     set(h1, 'Interpreter','none','fontsize',4)

    if nargin == 3
        printType = 'svg';
    end

    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize',[width height]);
    set(gcf,'PaperPosition',[0 0 width height]);
    box off;
    
    switch upper(printType)
        case 'PNG'
            print('-dpng',[fname '.png'])
        case 'PDF'
            print('-dpdf',[fname '.pdf'])
        case {'TIF', 'TIFF'}
            print('-dtiff','-r300',[fname '.tif'])
        case {'SVG'}
            print('-dsvg',[fname '.svg'])            
    end
    
    
    