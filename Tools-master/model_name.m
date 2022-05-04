function model = model_name(x)

    switch x 
        case 0
            model = 'OLS';
        case 1
            model = 'MNL';   
        case 2 
            model = 'Ordered Probit';        
        case 3 
            model = 'Poisson regression';            
        case 4 
            model = 'Negative Binomial regression';                
        case 5
            model = 'Zero Inflated Poisson regression';
        case 6
            model = 'Zero Inflated Negative Binomial regression';
    end
end