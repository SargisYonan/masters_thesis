function [ d, Xe, Ze, P ] = proximity( X, Xp, Z, variogram_model, var_params, range_extend )
    Ze = [];
    Xe = [];
    d = [];
    dists = [];
    P = [];
    
    for k = 1:length(Z)
        h = norm(X(k,:) - Xp);
        if (h == 0)
            h = 0.1;
        end
        
        if (h < var_params(1) + range_extend) % smaller than the range
            Ze = [Ze; Z(k)];
            Xe = [Xe; X(k,:)];
            d = [d; variogram_model(var_params, h)];
            dists = [dists; h];
        end
    end
    
    
    if isempty(Ze)
        for k = 1:length(Z)
            h = norm(X(k,:) - Xp);
            d = [d; variogram_model(var_params, h)];
        end
        
        Ze = Z;
        Xe = X;
    end
    
    for x = 1:length(Ze)
        for y = 1:length(Ze)
            h = norm(Xe(x,:) - Xe(y,:));
            
            if (h == 0)
                h = .1; % nugget
            end
            
            P(x,y) = variogram_model(var_params, h);
        end
    end
end

