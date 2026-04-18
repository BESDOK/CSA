function val = GD(PF, PFtrue)
% Generational Distance — yakınsama metriği (küçük = iyi)
% PF     : bulunan Pareto cephesi (K x M)
% PFtrue : gerçek Pareto cephesi referans noktaları (R x M)
K = size(PF, 1);
d = zeros(K, 1);
for i = 1:K
    diff = PFtrue - PF(i,:);             % R x M
    d(i) = min( sqrt(sum(diff.^2, 2)) ); % en yakın referans noktaya mesafe
end
val = sqrt(sum(d.^2)) / K;
end
