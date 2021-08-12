% Part of Problem set 4 - discrete probability theory
%{
Problem 3.3
1.-2. Checking theoretical values
%}
% Drawing 10000 values between 1,..9 
s = 9;
n = 10000;
r = randi([1 s],1,n);
trials = linspace(1,n,n);
[rmin,rmax] = bounds(r); % checking that numbers are within bounds
exp = mu(r,trials);
expt = mu_t(s,n);

% make plot
plot(trials,exp,'x')
hold on 
plot(trials,expt);

% ps: remember matlab has own mean and std functions
% checking if numerical mean and theoretical mean is the same
%checking if standard deviation is the same
sigma = std(r);
sigmat = stdevt(s);


mu_diff= abs(expt(end)-exp(end))/expt(end)
sigma_diff = abs(sigma - sigmat)/sigmat

%{
Comment: Number of trials to achieve 1% accuracy  (meaning error less than 1%)
standard deviation is more accurate than mean, so to reach an accurate mean
one need more trials (10^4).
%}

function sigma = stdevt(s) % theoritcal std
    sigma = 1/(2*sqrt(3))*sqrt(s^2 - 1);
end

function expvals = mu(r,s)
    expvals = 1./s.*cumsum(r);
end

function exp = mu_t(s,n)
exp = zeros(n,1) + 1/2*(s+1);
end