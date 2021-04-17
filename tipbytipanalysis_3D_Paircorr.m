

function [c12_3D_tips,set_comb,radii]=tipbytipanalysis_3D_Paircorr(SaveName,Tiploc,data1,data2,tipws,Resolution,pixelsize,randomtest)

rMax=tipws/2;
pw=tipws/2;

for t=1:size(Tiploc,2)

datarange(1)= Tiploc(1,t) - pw; % yloc
datarange(2)= Tiploc(1,t) + pw;
datarange(3)= Tiploc(2,t) - pw;
datarange(4)= Tiploc(2,t) + pw;


% z Tiploc (average among the data)
x1=data1.x;
y1=data1.y;

x2=data2.x;
y2=data2.y;

min_x=datarange(1);
max_x=datarange(2);
min_y=datarange(3);
max_y=datarange(4);


z1=data1.z;
z2=data2.z;

ind1=[];
ind1 =(x1>datarange(1) & x1<datarange(2)).* (y1>datarange(3) & y1<datarange(4));
z1_ws=z1(ind1==1);

ind2=[];
ind2 =(x2>datarange(1) & x2<datarange(2)).* (y2>datarange(3) & y2<datarange(4));
z2_ws=z2(ind2==1);

Tiploc(3,t)=mean(z1_ws);

N2_ws = sum(ind2);
x2_ws=x2(ind2==1);
y2_ws=y2(ind2==1);


datarange(5)= Tiploc(3,t) - pw;
datarange(6)= Tiploc(3,t) + pw;

min_z=datarange(5);
max_z=datarange(6);

ws=[(max_x-min_x),(max_y-min_y),(max_z-min_z)];


[c12_3D, radii,N1_ws,N2_ws]= pair_corr_cross_3D(data1, data2, datarange, Resolution, rMax, pixelsize);

c12_3D_tips(t,:)=c12_3D(1:end-1);


%% randomize test

if randomtest==1

   rand_data1.x= (max_x-min_x).*rand((N1_ws),1)+min_x;
   rand_data1.y= (max_y-min_y).*rand((N1_ws),1)+min_y; 
   rand_data1.z= (max_z-min_z).*rand((N1_ws),1)+min_z;


   rand_data2.x= (max_x-min_x).*rand((N2_ws),1)+min_x;
   rand_data2.y= (max_y-min_y).*rand((N2_ws),1)+min_y; 
   rand_data2.z= (max_z-min_z).*rand((N2_ws),1)+min_z;

   [c12_3D_Rand, radii,N1_ws,N2_ws]= pair_corr_cross_3D(rand_data1, rand_data2, datarange, Resolution, rMax, pixelsize);

   c12_3D_tips_rand(t,:)=c12_3D_Rand(1:end-1);

end

% 
end

Tiploc_3D=Tiploc;

set_comb.c12_tips= stat_combine_mat(c12_3D_tips,1,0,1);
set_comb.c12_tips_rand= stat_combine_mat(c12_3D_tips_rand,1,0,1);
set_comb.radii=radii;
set_comb.Tiploc_3D=Tiploc_3D;



figure
errorbar(radii(1:end-1),set_comb.c12_tips.mean,set_comb.c12_tips.ste);
hold on 
errorbar(radii(1:end-1),set_comb.c12_tips_rand.mean,set_comb.c12_tips_rand.ste);
legend ('data','random')
title ('pair-correlation')
xlabel('distance(nm)')
ylabel('c''')

saveas(gcf,['Paircorr_3D_' ,SaveName]);


   save(['Paircorr3D_',SaveName,'.mat'],'set_comb','-v7.3');

end