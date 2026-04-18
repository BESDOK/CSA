function val = IGD(PF, PFtrue)
% Inverted Generational Distance — yakınsama + yayılım (küçük = iyi)
% PF     : bulunan Pareto cephesi (K x M)
% PFtrue : gerçek Pareto cephesi referans noktaları (R x M)
R = size(PFtrue, 1);
d = zeros(R, 1);
for i = 1:R
    diff = PF - PFtrue(i,:);             % K x M
    d(i) = min( sqrt(sum(diff.^2, 2)) ); % gerçek noktadan bulunanlara en yakın mesafe
end
val = mean(d);
end
