% omitoption =1 % omit any set of data that has all zero 

function [result] = stat_combine_mat(MatData,dimention,omitoption,savematdata)
totalnum=size(MatData);

if dimention==1
    checkdimention=2;
elseif dimention==2
    checkdimention=1;
else
    error ('choose dimention between 1 or 2)')
end

% omit, remove zeros %
if omitoption==1
   omit_num=length(find(all(MatData==0 | isnan(MatData),checkdimention)));
  [removethis]=find(all(MatData==0 | isnan(MatData),checkdimention));
  if dimention==1
     MatData(removethis,:)=[];
  elseif dimention==2
     MatData(:,removethis)=[]; 
  end
else
omit_num=0;    
end


    result.totalnum=totalnum; 
    result.num=size(MatData,dimention);     
    result.omit_num=omit_num;
    result.dimention=dimention;
    result.mean=mean(MatData,dimention,'omitnan');
    result.std=std(MatData,[],dimention,'omitnan');
    result.ste=result.std./sqrt(result.num);
    
    result.median=median(MatData,dimention);
    result.quatile25=quantile(MatData,[0.25],dimention);
    result.quatile50=quantile(MatData,[0.5],dimention);
    result.quatile75=quantile(MatData,[0.75],dimention);

    if savematdata==1
     result.matdata=MatData;
    end

end



