%{

COLONY SEARCH ALGORITHM - MULTIOBJECTIVE ( CSA_MO )
Civicioglu,P.; Besdok, E.

Deb's rule (dominance-based selection)
kullanılarak multi-objective formuna dönüştürülmüş versiyon.

KULLANIM:
  out = CSA_MO(@zdt1, [], 100, 30,  0, 1, 500);
  out = CSA_MO(@zdt2, [], 100, 30,  0, 1, 500);
  out = CSA_MO(@zdt3, [], 100, 30,  0, 1, 500);
  out = CSA_MO(@zdt4, [], 100, 10, -5, 5, 500);
  out = CSA_MO(@zdt6, [], 100, 10,  0, 1, 500);

  plot(out.PF(:,1), out.PF(:,2), '.'); grid on




%{
================================================================================
  CSA_MO - ZDT1 TESTİ
================================================================================
  ZDT1: D=30, low=0, up=1
  Gerçek Pareto cephesi: f2 = 1 - sqrt(f1),  f1 in [0,1]  (sürekli, dışbükey)
================================================================================
%}

clear out out_CSA_MO

%% CSA_MO'yu ZDT1 üzerinde çalıştır
out = CSA_MO(@zdt1, [], 100, 30, 0, 1, 10000);

%% ZDT1'in gerçek Pareto cephesi: f2 = 1 - sqrt(f1)
f1t    = linspace(0, 1, 500)';
PFtrue = [f1t, 1 - sqrt(f1t)];

%% Görselleştirme
figure;
plot(out.PF(:,1), out.PF(:,2), 'r.', 'MarkerSize', 12); hold on
plot(PFtrue(:,1), PFtrue(:,2), 'k-', 'LineWidth', 1.2);
legend('bulunan','gerçek', 'Location', 'northeast');
grid on; xlabel('f_1'); ylabel('f_2');
title('CSA\_MO — ZDT1: bulunan vs. gerçek Pareto cephesi');

%% Metrikler
gd  = GD (out.PF, PFtrue);
igd = IGD(out.PF, PFtrue);
hv  = HV2D(out.PF, [1.1, 1.1]);

fprintf('\n=== CSA_MO - ZDT1 Sonuçları ===\n');
fprintf('GD  = %.6e\n', gd);
fprintf('IGD = %.6e\n', igd);
fprintf('HV  = %.6f\n', hv);
fprintf('|PF| = %d\n', size(out.PF,1));








%}

function out  =  CSA_MO(fnc , mydata , N , D , low , up ,Epk)
% Search Space limits
if numel(low)<D, low = low*ones(1,D); up = up*ones(1,D); end
% Initialization
T = 2  ;                      % Pattern matrix expanding factor or Colony (i.e., swarm)  size
initindex = 1:N;
moment = 0;
p0 = rand(T*N,D);     % Pattern matrix or Colony
% Initialize the Colony, p0
for i = 1:T*N
    for j = 1:D
        p0(i,j)  =  rand.*(up(j)-low(j))+low(j);
    end
end
fitp0  =  feval(fnc,p0,mydata);       % (T*N) x M  (M = amaç sayısı)
% Iterative Search Phase
for epoch  =  1:Epk
    % Selection of predators, p
    while 1
        index = randperm(T*N,N); % uniform selection
        if sum(index == initindex) == 0, initindex = index; break; end
    end
    p = p0(index,:);
    fitp = fitp0(index,:);             % N x M
    % Direction scaling factor, scale
    if rand < rand, c  =  1; else, c  =  D; end
    if rand<rand
        scale  =  ( rand(N,c)-0.50)  ./ (rand(N,c)-0.50);
    else
        t = [-1 1];  scale = sign(rand(N,1)-0.50) .* levy_dist(randi([2 5],N,1),randi(10,N,1).^t(randi(2)));
    end
    % Morphogenesis (i.e., mutation + crossover) control matrix, m
    m = zeros(N,D);
    for j = 1:N
        ind = randperm(D);
        k = abs( randi([0 1]) - rand^randi([2 10]) );
        b = ind( 1: ceil(k*D) );
        m(j,b) = 1;
    end
    % Setting of evolutionary direction vector, dx
    while 1, v1  =  randperm(N); v2  =  randperm(N);  if sum(v1  ==  (1:N)) ==  0 && sum(v2  ==  v1) ==  0 && sum(v2  ==  (1:N)) ==  0, break;  end, end
    [~,index0] = sort(fitp(:,1));      % MO: ilk amaca göre sıralama
    v = randi(3); % Random selection of interaction model
    switch v
        case 1, dx  =  p(v2,:) - p(v1,:)    ;                                % Bilateral-Bijective Interaction
        case 2, dx  = p(v1,:) - p   ;                                          % Bijective Interaction
        case 3, dx  = p( index0(randi([1 ceil(N/5) ])),:) - p  ;   % Swarmmic Interaction
    end
    % morphogenesis pattern matrix, px
    s = (rand(N,1)-0.50) .* rand(N,1).^randi( [2 10],1) ;
    px = p + scale .* m.* dx  + s .* moment ;
    % boundary control
    for k = 1:N
        for l = 1:D
            if px(k,l)<low(l), px(k,l) = low(l) + rand.^randi([ 1 5 ] ,1)*( up(l) - low(l) ); end
            if px(k,l)>up(l),  px(k,l) = up(l)  + rand.^randi([ 1 5 ],1)*( low(l)-  up(l) ); end
        end
    end
    % Deb's rule (dominance-based selection)
    fitpx  =  feval(fnc,px,mydata);
    ind = all(fitpx <= fitp, 2) & any(fitpx < fitp, 2);
    % Update Clan
    p(ind,:)  =  px(ind,:);
    fitp(ind,:)  =  fitpx(ind,:);
    % Update Colony
    p0(index,:)  =  p;
    fitp0(index,:)  =  fitp;

    %% Pareto cephesi güncelle
    ndMask   = paretoFront(fitp0);
    PF       = fitp0(ndMask,:);
    PS       = p0(ndMask,:);

    out.PF = PF;   % Pareto Front (amaç uzayı)
    out.PS = PS;   % Pareto Set   (karar uzayı)
    assignin('base','out',out)
    fprintf('%s | %5.0f --> |PF|=%4d\n', func2str(fnc), epoch, size(PF,1));
    % Update momentum
    moment =  (abs(randi([0 1],N,1))-m).*dx;

end
end

function nd = paretoFront(F)
% Pareto-optimal (non-dominated) kümesinin mantıksal maskesi
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

function x  =  levy_dist(alpha, beta)
% Generate a random number from a Lévy distribution with characteristic
% exponent `alpha` and scale parameter `beta`.

% First, generate a standard Cauchy random variable.
z  =  rand() + 1;

% Then, generate a random number from the Gamma distribution with shape
% parameter `alpha` and scale parameter 2.
w  =  gamrnd(alpha, randi([2 5]));

% Finally, compute the Lévy random variable `x`.
x  =  beta * z ./ w.^(1 ./ alpha);
end