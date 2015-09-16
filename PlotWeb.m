function PlotWeb(matrix)
    PosX                = []; %X coordinate of each cell in a unit square
    PosY                = []; %Y coordinate of each cell in a unit square
    nRows               = []; %row number
    nCols               = []; %col number
    CellColor           = [1 1 1]; %Cell color
    BackColor           = [0 0 128]/255; %Back color
    Margin              = 0.12; %Margin between cels
    Axis                = 'equal'; %Axis
    BorderColor         = 'none';
    BorderLineWidth     = 0.0005;

    [nRows nCols] = size(matrix);
    PosX = zeros(nRows,nCols);
    PosY = zeros(nRows,nCols);
    
    PosX = repmat(1:nCols, nRows, 1);
    PosY = repmat(((nRows+1)-(1:nRows))',1,nCols);
    
    DrawBack(BackColor); %Draw the back color

    for i = 1:nRows
        for j = 1:nCols
            if(matrix(i,j) > 0)
                DrawCell(i,j,CellColor);
%                 if(i==1)
%                     DrawCell(i,j,[0 1 0]);
%                 end
            end
        end
    end

    ApplyBasicFormat;
    
    function DrawBack(color)
        rec = rectangle('Position',[0.5-Margin-0.01,0.5-Margin, nCols+2*Margin,nRows+2*Margin]);
        axis(Axis);
        xlim([0.5-1.1*Margin,0.5+nCols+Margin]);
        ylim([0.5-Margin,0.5+nRows+Margin]);
        
        set(gca,'box','on');
        set(gcf,'DefaultAxesLineWidth',0.0005);
        set(rec,'FaceColor',color);
        set(rec,'EdgeColor',BorderColor);
        set(rec,'LineWidth',BorderLineWidth);
    end

    function DrawCell(i,j,color)

        rec = rectangle('Position',[PosX(i,j)-0.5+Margin,PosY(i,j)-0.5+Margin,1-2*Margin,1-2*Margin]);
        set(rec,'FaceColor',color);
        set(rec,'EdgeColor','none');

    end

    function ApplyBasicFormat
       
        set(gca,'xticklabel',[]);
        set(gca,'yticklabel',[]);
        set(gca,'YTick',[]);
        set(gca,'XTick',[]);

    end
end