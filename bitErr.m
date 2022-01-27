function [num_err, err_rate] = bitErr(x,y)
    % compare x and y to get # of symbol errors (bit errors for binary data)
    num_err = length(find(x ~= y));
    err_rate= vpa(num_err/size(x,1)/size(x,2));
end