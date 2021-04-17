% colocalization analysis between image A and B
% R: Pearson’s correlation coefficient (R’), overlap coefficient (MOC), and Manders’ correlation coefficients 1 and 2 (M1, and M2)
% Reference: Manders, E. M. M., Verbeek, F. J. & Aten, J. A. Journal of Microscopy 169, 375-382, doi:10.1111/j.1365-2818.1993.tb03313.x (1993).
% test with loading the 'colocalization_example_ch1.mat','colocalization_example_ch2.mat' file
% Thresh_Option 0: no threshold ; 1: set a threshold using multithresh

%% loading example data 

%load('colocalization_example_ch1.mat')
%load('colocalization_example_ch2.mat')
%imshowpair(ch2,ch1)

function [CorrResult]=colocalization_analysis(ch1,ch2,Thresh_Option)


if Thresh_Option==0
    [R,MOC,k1,k2,rsqr,M1,M2] = colocalization(ch1,ch2);
    
elseif Thresh_Option==1
    
    [ThreshOut_ch1,~,ThreshImage_ch1] = threshold_image(ch1);
    [ThreshOut_ch2,~,ThreshImage_ch2] = threshold_image(ch2);
    [R,MOC,k1,k2,rsqr,M1,M2] = colocalization(ThreshImage_ch1,ThreshImage_ch2);
else
    error('choose 0 or 1')
end

%CorrResult= [R,MOC,k1,k2,rsqr,M1,M2];
CorrResult.R=R;
CorrResult.MOC=MOC;
CorrResult.k1=k1;
CorrResult.k2=k2;
CorrResult.rsqr=rsqr;
CorrResult.M1=M1;
CorrResult.M2=M2;
end

function [Thresh,BinaryImage,ThreshImage] = threshold_image(image)

level = multithresh((image),2);
Thresh=level(1);

ThreshImage=image-Thresh;
ThreshImage(ThreshImage<0)=0;
ThreshImage=round(ThreshImage,0);
BinaryImage= imbinarize(image,Thresh);
if max(BinaryImage(:))==min(BinaryImage(:))
    BinaryImage=ones(size(image));
end
end

function[R,MOC,k1,k2,rsqr,M1,M2] = colocalization(A,B)

Ai=A(B>0);
Bi=B(A>0);

BW1=(A>0);
BW2=(B>0);
BW3=BW1+BW2;

A(BW3==0)=nan;
B(BW3==0)=nan;

AA=A.*A;
BB=B.*B;
AB=A.*B;
mA=mean(A(:),'omitnan');
mB=mean(B(:),'omitnan');

%% MOC: Manders Overall Overlab Coefficienct
a=(A-mA);
b=(B-mB);
ab=a.*b;
aa=a.*a;
bb=b.*b;

R=sum(ab(:),'omitnan')/(sqrt(sum(aa(:),'omitnan'))*sqrt(sum(bb(:),'omitnan'))); %% pixel by pixel S (Colocalization)
MOC=sum(AB(:),'omitnan')/sqrt((sum(AA(:),'omitnan')).*(sum(BB(:),'omitnan'))); %% pixel by pixel S (Colocalization)

k1= sum(AB(:),'omitnan')/sum((AA(:)),'omitnan');
k2= sum(AB(:),'omitnan')/sum((BB(:)),'omitnan');
rsqr=k1*k2;

M1= sum(Ai(:),'omitnan')/sum(A(:),'omitnan');
M2= sum(Bi(:),'omitnan')/sum(B(:),'omitnan');
end