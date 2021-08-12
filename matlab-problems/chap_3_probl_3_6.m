% Discrete probability theory: problem 3.6
% generalized dice

%}
% Drawing 10000 values between 1,..9 
s = 2;  % how many sides on die has
N = 200; % number of die
n = 10000; % number of experiments

sum_ = zeros(1,n); % important that dimensions the same
for i = 1:n
    r = randi([1 s],1,N);
    sum_(i) = sum(r);
end

sum_max = 9*N;
% obtaining numerical mean, variance and standard deviation
trials = linspace(1,n,n);
exp = mu(sum_,trials);

% make plot
plot(trials,exp,'x');


%sigma and mean
var_ = var(sum_)
sigma_t = stdevt(s,N)
sigma = std(sum_)
mu_num = exp(end)
mu_th = mu_t(s,N) 

% histogram (from uniform dice to non-uniform sum :))
hist(sum_)

% check theoretical if time
mu_diff= abs(mu_num-mu_th)/mu_th
sigma_diff = abs(sigma - sigma_t)/sigma_t

%{
Comment: Note that relative standard deviation sigma_x/mu_x goes as
1/swrt(n), where n is number of sample data/experiments

Results:
a) two dice, 10 sides
n = 1000
var_ = 16.3143
sigma =  4.0391
mu_num = 10.9450
Note: variance increases with a larger samplesize, but we get closer to
real mean

b) 10 dice, 20 sides



%}


function sigma = stdevt(s,N) % theoritcal std
    sigma = 1/(2*sqrt(3))*sqrt(s^2 - 1)*sqrt(N);
end


function exp = mu_t(s,N)
exp =  1/2*(s+1)*N;
end

function expvals = mu(r,s)
    expvals = 1./s.*cumsum(r);
end


