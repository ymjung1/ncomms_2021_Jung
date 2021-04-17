function [SR_Seg]=calculate_SR(ch1,ch2,BW,ch1_part,ch2_part)


ch1_seg=ch1(BW==1);
ch2_seg=ch2(BW==1);

meanch1=mean(ch1_seg);
meanch2=mean(ch2_seg);
mean_ch12=[meanch1 meanch2];



s_ch1_hist=sqrt((sum((ch1_part-meanch1).^2,1))/length(ch1_part));
s_ch2_hist=sqrt((sum((ch2_part-meanch2).^2,1))/length(ch2_part));


Xhist=[ch1_part, ch2_part] ; %(part of X data)
xc = Xhist - mean_ch12;  % subtract mean

c_hist = (xc' * xc) ./ size(Xhist,1);

coeff_hist=c_hist/(s_ch1_hist*s_ch2_hist);
SR_Seg=coeff_hist(2);


end