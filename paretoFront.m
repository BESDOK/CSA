function nd = paretoFront(F)
% Pareto-optimal (non-dominated) kümesinin mantıksal maskesi
% F  : K x M matris (K nokta, M amaç)
% nd : K x 1 mantıksal vektör — true ise satır non-dominated
n = size(F,1);
nd = true(n,1);
for i = 1:n
    if nd(i)
        for j = 1:n
            if i ~= j && all(F(j,:) <= F(i,:)) && any(F(j,:) < F(i,:))
                nd(i) = false;
                break;
            end
        end
    end
end
end