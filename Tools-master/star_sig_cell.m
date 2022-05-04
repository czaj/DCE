function star_array = star_sig_cell(P)

star_array = cell(size(P));
    
for j = 1:size(P,2)
    for i = 1:size(P,1)
        if P(i,j) <= 0.01
            star_array(i,j) = {'***'};
        elseif P(i,j) <= 0.05
            star_array(i,j) = {'** '};
        elseif P(i,j) <= 0.1
            star_array(i,j) = {'*  '};
        else
            star_array(i,j) = {'   '};
        end
    end
end
