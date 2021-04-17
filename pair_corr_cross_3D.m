% dr: increment for increasing radius 
% rMax: largest radius

function [g, radii,N1_ws,N2_ws]= pair_corr_cross_3D(data1,data2, myrange, dr, rMax, pixelsize)

x1=data1.x;
y1=data1.y;
z1=data1.z;

x2=data2.x;
y2=data2.y;
z2=data2.z;


%%   RECT=[];

min_x=myrange(1);
max_x=myrange(2);
min_y=myrange(3);
max_y=myrange(4);
min_z=myrange(5);
max_z=myrange(6);

ws=[(max_x-min_x),(max_y-min_y),(max_z-min_z)];

dr =dr/pixelsize;
rad_bins = 0:dr:rMax; 
in_radii = rad_bins(1:end-1); 
out_radii = rad_bins(2:end); 
mean_r = (in_radii+out_radii)/2; 
radii=mean_r*pixelsize;
volumns = (4 / 3)*pi*(out_radii.^3 - in_radii.^3); 

ind1 =(x1>myrange(1) & x1<myrange(2)).* (y1>myrange(3) & y1<myrange(4)).* (z1>myrange(5) & z1<myrange(6));

N1_ws = sum(ind1);
x1_ws=x1(ind1==1);
y1_ws=y1(ind1==1);
z1_ws=z1(ind1==1);


ind2 =(x2>myrange(1) & x2<myrange(2)).* (y2>myrange(3) & y2<myrange(4)).* (z2>myrange(5) & z2<myrange(6));

N2_ws = sum(ind2);
x2_ws=x2(ind2==1);
y2_ws=y2(ind2==1);
z2_ws=z2(ind2==1);


 if N1_ws < 1 
     error( " particle must be more than 1 ")
 else
     
     
    g = zeros(N2_ws, length(rad_bins)-1);   
    totalDensity = (N2_ws)/(ws(1)*ws(2)*ws(3));

         cnt=0;
     for p=1:N1_ws
         distance = sqrt((x1_ws(p) - x2_ws).^2 + (y1_ws(p) - y2_ws).^2 + (z1_ws(p) - z2_ws).^2);
         [result] = (hist(distance,rad_bins(1:end-1))); 
         g(p,:) = result/totalDensity;
     end
    

    g = mean(g,1)./volumns;

 end
 end 
