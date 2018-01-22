function [ g ] = kernel( h, range, sill, kernel_model )
% computes the statistical kernel passed in on the values passed in

    function g = gaussian_kernel(h, a, sill)
        g = sill * ones(size(h)) - sill * exp(-3 .* h.^2 ./ a^2);
    end

    function g = spherical_kernel(h, a, sill)
        g = zeros(size(h));
        for i = 1:max(size(g))
            if h(i) <= a
                g(i) = sill * 1.5 * (h(i) / a) - 0.5 * sill * (h(i) / a).^3;
            else
                g(i) = sill;
            end
        end
    end

    function g = exponential_kernel(h, a, sill)
        g = sill * ones(size(h)) - sill * exp(-3 .* h.^2 ./ a^2);
    end

    switch(kernel_model)
        case 'gaussian'
            g = gaussian_kernel(h, range, sill);
        case 'normal'
            g = gaussian_kernel(h, range, sill);
        case 'exp'
            g = exponential_kernel(h, range, sill);
        case 'exponential'
            g = exponential_kernel(h, range, sill);
        case 'sph'
            g = spherical_kernel(h, range, sill);
        case 'spherical'
            g = spherical_kernel(h, range, sill);
        otherwise
            g = NaN;
            disp('Error: pass a valid kernel option into kernel_model')
    end
    
end


