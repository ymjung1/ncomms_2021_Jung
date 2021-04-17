% draw scalebar (nm scale) starting from 10% of range x and range y
% eg: addscale5_Auto_range(1000,[1 1 1],3)
% Exfactor,markerlocation, positionfactor=15 

function[ret]=add_scalebar(markersize,scalenarcolor,scalethinkness,pixelsizenm,Exfactor,markerlocation)

markersize=markersize*Exfactor;
range_x= get(gca,'xlim');
range_y= get(gca,'ylim');
ROIPix_Y=range_y(2)-range_y(1)+1;
ROIPix_X=range_x(2)-range_x(1)+1;
%starty=(ROIPix_Y/10)+range_y(1); 
scalerange=markersize/pixelsizenm;
Markerx=linspace(0,scalerange);
Markery=ones(size(Markerx));
positionfactor=15;
switch markerlocation
    case 1 % 'northwest'
starty=(ROIPix_Y/positionfactor)+range_y(1); %%  fixed for image axis
startx=(ROIPix_X/positionfactor)+range_x(1);
    case 2 %'southwest'
starty=(ROIPix_Y/positionfactor)+range_y(1); %%  fixed for image axis
startx=ROIPix_X-(ROIPix_X/positionfactor)-range_x(1)-scalerange;
    case 3 %'northeast'
starty=ROIPix_Y-((ROIPix_Y/positionfactor)+range_y(1)); %%  fixed for image axis
startx=(ROIPix_X/positionfactor)+range_x(1);
    case 4 %'northwest'
starty=ROIPix_Y-((ROIPix_Y/positionfactor)+range_y(1)); %%  fixed for image axis
startx=ROIPix_X-(ROIPix_X/positionfactor)-range_x(1)-scalerange;
end


Markerx=startx+Markerx;
Markery=starty*Markery;
scalebar=line(Markerx,Markery,'Color',scalenarcolor,'LineWidth',scalethinkness); 


end

