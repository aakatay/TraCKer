function y = linearfit(x,v)
% similiar to poloyfit but a linear betw, data points

    j = 1;
    for i = 1:numel(x)-1
        yix = x(i):x(i+1);
        ix = yix-x(i);
        y(yix) = v(i)+ (v(i+1)-v(i))/(x(i+1)-x(i))*ix;
    end
end