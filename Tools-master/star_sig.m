function star_array = star_sig(P)

star_array = [];
for i = 1:length(P)
    if P(i) <= 0.01 
        star_array = [star_array; '***   '];
    elseif P(i) <= 0.05 
        star_array = [star_array; '**    '];
	elseif P(i) <= 0.1 
        star_array = [star_array; '*     '];
    else 
        star_array = [star_array; blanks(6)];
    end
end
