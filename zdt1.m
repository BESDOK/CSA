function f = zdt1(P, userdata)
% ZDT1 - convex Pareto front
% Önerilen: D = 30, x in [0,1]
[N, D] = size(P);
f1 = P(:,1);
g  = 1 + 9 * sum(P(:,2:D), 2) / (D - 1);
h  = 1 - sqrt(f1 ./ g);
f2 = g .* h;
f  = [f1, f2];
end
