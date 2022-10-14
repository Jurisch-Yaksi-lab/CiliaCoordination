function [resultData]=jnh_FreqFilter(data,SDthreshold,Fs)
%Copyright 2017, JN Hansen
%data must be a two dimensional matrix
%% filter by SD value - TODO only filter out if out of range
resultData = zeros(size(data,1), size(data,2));
for i=1:size(data,1)
    for j=1:size(data,2)
        if(j==1||i==1||j==size(data,2)||i==size(data,1))
        elseif(std2(data(i-1:i+1,j-1:j+1))<SDthreshold)
            resultData(i-1:i+1,j-1:j+1)=ones(3,3);
        end
    end
end

for i=1:size(data,1)
    for j=1:size(data,2)
        if(resultData(i,j)~=0 && data(i,j)<Fs/2)
            resultData(i,j) = data(i,j);
        else
            resultData(i,j) = NaN;
        end
    end
end
