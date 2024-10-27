function errtab = summarize_errors(errcell, print)

N_runs = length(errcell);
NT_sim = size(errcell{1},1);
N_eqn = size(errcell{1},2);

errmat = zeros( N_runs * NT_sim, N_eqn );

for ii=1:N_runs
	errmat( (ii-1)*NT_sim + (1:NT_sim), : ) = errcell{ii};
end

avg_err=mean(abs(errmat))';
med_err=median(abs(errmat))';
p75_err=prctile(abs(errmat),75)';
p95_err=prctile(abs(errmat),95)';
p99_err=prctile(abs(errmat),99)';
p995_err=prctile(abs(errmat),99.5)';
max_err=max(abs(errmat))';
errtab=table(avg_err,med_err,p75_err,p95_err,p99_err,p995_err,max_err);
errarr=table2array(errtab);

if print
	disp(' ');
	disp('-----------------------------------------------');
	disp('Average and maximum Euler equation error');
	fprintf('Equ.no.\t\tAvg.\t\tMed.\t\tp75\t\t\tp95\t\t\tp99\t\t\tp99.5\t\tMax.\n');
	for s=1:length(avg_err)
		fprintf('%d\t\t\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',s,errarr(s,1),errarr(s,2),errarr(s,3), ...
			errarr(s,4),errarr(s,5),errarr(s,6),errarr(s,7));
	end
end