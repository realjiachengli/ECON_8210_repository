function exst = mcmc(T,mtrans,permind,start,shmat)

if nargin<5
	shmat = lhsdesign(T,1,'criterion','correlation');
end

exst = zeros(T,1);
exst(1) = start;

for t=1:T-1
    % next period's exog. state
	exstperm=permind(exst(t));
    transprob=cumsum(mtrans(exstperm,:));
    nextst=find(transprob-shmat(t)>0,1,'first');
    exst(t+1)=nextst;
	end
end