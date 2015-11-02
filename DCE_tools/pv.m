function result = pv(m,s)

if any(s < 0) || ~isreal(s)
    result = zeros(size(m));
    for i = 1:length(m)
        if s(i) < 0 || ~isreal(s(i))
            result(i) = NaN;
        else
            result(i) = (1-normcdf(abs(m(i))./real(s(i)),0,1))*2;
        end
    end
else
    result = (1-normcdf(abs(m)./s,0,1))*2;
end
