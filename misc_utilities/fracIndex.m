function fracIndx = fracIndex(array,x)
% Finds the fractional index of 'x' in the array 'array'. 
% For example, if x = 4.5, and array = [4, 5, 6, 7], then fracIndex = 1.5.
% x can be an array also, in which case fracIndex is an array the same size
% of x. 
% Using linear interpolation between nearest points.
fracIndx = zeros(1,length(x));
 for idx = 1:length(x)
    if x >= array(end)
        fracIndx(idx) = length(array);
    elseif x(idx) <= array(1)
        fracIndx(idx) = 1;
    else
        a = find(array <= x(idx));
        a = a(end);
        b = find(array > x(idx));
        b = b(1);
        fracIndx(idx) = a+(x(idx)-array(a))/(array(b)-array(a));
    end
    
end