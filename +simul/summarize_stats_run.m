function [statsout_run, freq_exog, addstats] = summarize_stats_run( simseries, ~ )

[~,dispnames,nvars] = simul.get_simseries_array(simseries);
gdp_idx = find(dispnames=="Y"); 
zp_idx = find(dispnames=="zp");

% first for all periods, then separately for good and bad states
smpsel_exog{1} = @(simseries) true(size(simseries,1),1);
smpsel_exog{2} = @(simseries) ~simseries.isProfligacy & ~simseries.isAusterity;
smpsel_exog{3} = @(simseries) simseries.isProfligacy;
smpsel_exog{4} = @(simseries) simseries.isAusterity;
N_subsamples = length(smpsel_exog);

statsout_run=cell(N_subsamples,1);
addstats=cell(N_subsamples,1);

freq_exog = zeros(N_subsamples,1);

[~, colnames] = simul.estimate_sample_moments(zeros(0,nvars),[zp_idx,gdp_idx]);

for j=1:N_subsamples 
	smpsel = smpsel_exog{j}(simseries);
	% convert single columns to a matrix
	simarray = simul.get_simseries_array(simseries);

	% subsample
	simtmp = simarray(smpsel,:);
	freq_exog(j)=size(simtmp,1);
	
	statstmp = simul.estimate_sample_moments(simtmp,[zp_idx,gdp_idx]);
	
	statsout_run{j}=array2table(statstmp,'RowNames',dispnames,'VariableNames',colnames);
	addstats{j} = computeAdditionalStats( simseries, smpsel );
end