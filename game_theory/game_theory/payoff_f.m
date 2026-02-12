function [w1, w2] = payoff(s1, s2, n, U)
w1 = 0;
w2 = 0;
for i=1:n
    w1 = w1 + U(s1(i), s2(i));
    w2 = w2 + U(s2(i), s1(i));
    if s1(i)==2 || s2(i)==2
        break
    end
end
