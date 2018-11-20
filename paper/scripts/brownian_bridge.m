% X0 = [ceil(100 * rand(2,1)); 0];
X0 = [50, 6, 0]';

mck = 15;
N = 2;
KN = [[5 80];
     [95 90]];

% KN = ceil(100 * rand(N, 2));

alpha = 5;
varw = alpha/2;

figure();

for i = 1 : N
    px = KN(i,1);
    py = KN(i,2);

    plot([X0(1) px], [X0(2) py], '-r', 'LineWidth', 2)
    hold on;
    
    for m = 1 : mck
    T = [];
    T = X0;
    T(3) = atan2d(py-X0(2), px-X0(1));
    
    D = norm([px py]' - X0(1:2));
    for k = 2: 2*ceil(D/alpha)
        thetak = T(3,k-1);

        Dx = abs(px - T(1, k-1));
        varx = abs((cosd(thetak) * Dx) / (cosd(thetak) + Dx));

        Dy = abs(py - T(2, k-1));
        vary = abs((sind(thetak) * Dy) / (sind(thetak) + Dy));

        xk = T(1,k-1) + alpha * cosd(thetak) + sqrt(varw)*randn;
        yk = T(2,k-1) + alpha * sind(thetak) + sqrt(varw)*randn;
        thetak = atan2d(py-yk, px-xk);

        T(:,k) = [xk, yk, thetak]';

        if (norm([xk yk] - [px py]) <= alpha)
             T(:,k+1) = [px, py, 0]';
            break;
        end
    end

    plot(T(1,:), T(2,:), 'k', 'LineWidth', .8)
    hold on;
    plot(px, py, 'ko', 'LineWidth', 50);
    hold on;
   
    size(T)
    end

end

plot(X0(1), X0(2), 'go', 'LineWidth', 50)
hold off;
axis([1 100 1 100])

