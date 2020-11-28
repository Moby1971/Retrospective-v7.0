function window = respiratory_window(rlocs,mean_resp,resp_percentage,tr,nr_klines)

% this function creates an array (time line) of rectangular boxes of 0's and 1's around the detected respiratory signals
% later on 1 means that there is a respiration, for which the k-lines will be discarded

resp_a = round(0.5*(resp_percentage/100)*(60/mean_resp)*1000/tr);   % 1/2 window width around respiratory peak locations

temp_window=zeros(nr_klines,1);   % fill array with zeros

for i=1:size(rlocs,2)
    center = round(rlocs(1,i));
    for j=center-resp_a:center+resp_a
        if (j>0) && (j<= nr_klines)
            temp_window(j)=1;
        end
    end
end

window = temp_window;

end