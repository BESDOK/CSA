function val = HV2D(PF, refPoint)
% Hypervolume — 2 amaçlı problem için (büyük = iyi)
% PF       : bulunan Pareto cephesi (K x 2)
% refPoint : referans (nadir) nokta, örn ZDT1/2/3/4/6 için [1.1 1.1] uygundur
%
% Referans noktayı domine etmeyen satırlar atılır; sonra sol->sağ sıralanıp
% dikdörtgen toplamıyla alan hesaplanır.
mask = all(PF < refPoint, 2);
F = PF(mask, :);
if isempty(F), val = 0; return; end
F = sortrows(F, 1);          % f1'e göre artan sırala
val = 0;
prevF2 = refPoint(2);
for i = 1:size(F,1)
    w = refPoint(1) - F(i,1);
    h = prevF2 - F(i,2);
    if h > 0
        val = val + w * h;
        prevF2 = F(i,2);
    end
end
end
