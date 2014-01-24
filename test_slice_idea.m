function samples=test_slice_idea(nsamples)

contprob=@(x) normpdf(x); 
theta=1;
w=3.0;
% f(x) = contprob(x) + w delta_theta(x)
delta=1;

samples=slicesample(0.0, nsamples, 'pdf', @(x)effective_f(x,contprob,theta,w,delta));

function p=effective_f(x,contprob,theta,w,delta)
    % assume theta>0
    if x < theta
        p=contprob(x); 
    elseif x < (theta+delta)
        p=w/delta;
    else
        p=contprob(x-delta); 
    end
end

end