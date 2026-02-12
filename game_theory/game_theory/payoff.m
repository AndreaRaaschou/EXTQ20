function [w1, w2] = payoff(p1, p2, n, U)
% set empty start-values    
w1 = 0;
w2 = 0;

% for the amount of strategies in the sequence
for i=1:n
    % add points from U to w1 and w2
    w1 = w1 + U(p1(i), p2(i));
    w2 = w2 + U(p2(i), p1(i));
    if p1(i) == 2 || p2(i) == 2
        break; % Exit the loop if either player uses strategi 2
    end
end
end

