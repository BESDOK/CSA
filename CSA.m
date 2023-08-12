%{


COLONY SEARCH ALGORITHM  ( CSA )
Civicioglu,P.; Besdok, E. 

12.August.2023

%}

function out  =  CSA(fnc , mydata , N , D , low , up ,Epk)
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
fitp0  =  feval(fnc,p0,mydata);
% Iterative Search Phase
for epoch  =  1:Epk
    % Selection of predators, p
    while 1
        index = randperm(T*N,N); % uniform selection
        if sum(index == initindex) == 0, initindex = index; break; end
    end
    p = p0(index,:);
    fitp = fitp0(index);
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
    [~,index0] = sort(fitp);
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
    % Greedy Selection
    fitpx  =  feval(fnc,px,mydata);
    ind  =  fitpx < fitp;
    % Update Clan
    p(ind,:)  =  px(ind,:);
    fitp(ind)  =  fitpx(ind);
    % Update Colony
    p0(index,:)  =  p;
    fitp0(index)  =  fitp;
    % Update global solutions
    [gmin,indbest]  =  min(fitp0);
    gbest  =  p0(indbest,:);
    % Report and return solutions
    [out.gmin , out.gbest]  =  deal(gmin , gbest  );
    assignin('base','out',out)
    fprintf('%3.0f-> %5.16f \n',epoch,gmin)
    % Update momentum
    moment =  (abs(randi([0 1],N,1))-m).*dx;
    % stopping
    if gmin<1e-12, return, end

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



