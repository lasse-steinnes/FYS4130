% Sums of 0s and 1s : modify program for binomial distribution
% Discrete probability theory: problem 3.6
% generalized dice

%}
% Drawing either value of 1 or 0 and taking the sum 
N = 10; % number of particles 
n = 1000; % number of experiments
p = 0.5;

sum_ = zeros(1,n); % important that dimensions the same
for i = 1:n
    r = binornd(N,p);
    sum_(i) = sum(r);
end

% obtaining numerical mean, variance and standard deviation
trials = linspace(1,n,n);
exp = mu(sum_,trials);


%sigma and mean
%var_ = var(sum_)
sigma = std(sum_)
sigma_t = std_t(p,N)
mu_num = exp(end)
mu_th = mu_t(p,N)

% plot
%fit distribution to data
sum_range = linspace(0,N,n);
gauss_b = gbinom(sum_range,mu_th,sigma_t)*n; % of total n tries, here is the distribution of sums
% histogram (binomial plot))
hist(sum_)
hold on
plot(sum_range,gauss_b, "r-")

% check theoretical if time
sigma_diff = abs(sigma_t-sigma)/sigma_t
mu_diff = abs(mu_num - mu_th)/mu_th

%{
Comment: Note the distinction between 
probability distribution (P(n|N))and probability (p) for each outcome!
%}

function expval = mu_t(p,N) % theoretical for binomial distribution
    expval = p*N;           % note, we need several "experiments/trials" 
end                         % to get this numerically
                            % exp and std for sum of binomial identically
                            % distributed integers (1/0).

function sigma = std_t(p,N)
    sigma = sqrt(p*(1-p)*N);
end 


function expvals = mu(r,s)
    expvals = 1./s.*cumsum(r);
end

function g = gbinom(x,mu,sigma)
g = 1./sqrt(2.*sigma.*pi.^2).*exp(-((x-mu).^2)./(2*sigma.^2));
end
