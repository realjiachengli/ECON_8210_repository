function [statsout_exog, stderrs_exog, dispnames, freq_means, addtables] = summarize_stats( simseries, ~ )

N_runs = max(simseries.run);
NT_sim = max(simseries.per);
[~,dispnames,nvars] = simul.get_simseries_array(simseries);
gdp_idx = find(dispnames=="Y"); 
zp_idx = find(dispnames=="zp");

% first for all periods, then separately for good and bad states
smpsel_exog{1} = @(simseries) true(size(simseries,1),1);
smpsel_exog{2} = @(simseries) ~simseries.isProfligacy & ~simseries.isAusterity;
smpsel_exog{3} = @(simseries) simseries.isProfligacy;
smpsel_exog{4} = @(simseries) simseries.isAusterity;
N_subsamples = length(smpsel_exog);

statsout_exog=cell(N_subsamples,1);
stderrs_exog=statsout_exog;

freq_exog = zeros(N_subsamples,N_runs);

[~, colnames] = simul.estimate_sample_moments(zeros(0,nvars),[zp_idx,gdp_idx]);

for j=1:N_subsamples 
	smpsel = smpsel_exog{j}(simseries);
	% convert single columns to a matrix
	simarray = simul.get_simseries_array(simseries);

	statstmp=zeros(nvars,length(colnames),N_runs);
	simcell = arrayfun(@(ii)simarray( smpsel & simseries.run==ii, : ),1:N_runs, ...
		'UniformOutput', false);
	%for ii=1:N_runs
	parfor ii=1:N_runs
		% subsample
		simtmp = simcell{ii};
		freq_exog(j,ii)=size(simtmp,1);
		
		slice = simul.estimate_sample_moments(simtmp,[zp_idx,gdp_idx]);
		statstmp(:,:,ii) = slice;
	end
	
	tmpmeans=mean(statstmp,3,'omitnan');
	statsout_exog{j}=array2table(tmpmeans,'RowNames',dispnames,'VariableNames',colnames);
	
	tmpstderrs=std(statstmp,[],3,'omitnan') / sqrt(N_runs); 
	stderrs_exog{j}=array2table(tmpstderrs,'RowNames',dispnames,'VariableNames',colnames);
end

% Additional variables
addcell = cell(N_runs, N_subsamples);
simcell = arrayfun(@(i)simseries(simseries.run==i,:),1:N_runs,'UniformOutput',false);
%for ii=1:N_runs
parfor ii=1:N_runs
	row = cell(1,N_subsamples);
	simrun = simcell{ii};
	for j=1:N_subsamples
		smpsel = smpsel_exog{j}(simrun);
		row{j} = computeAdditionalStats( simrun, smpsel );
	end
	addcell(ii,:)=row;
end

addtables = cell(N_subsamples,1);
for j=1:N_subsamples
	% Deal with additional stats
	tmpadd = vertcat( addcell{:,j} );
	tmpaddmeans = varfun( @(x)nanmean(x), tmpadd );
	tmpaddstderrs = varfun( @(x)nanstd(x)/sqrt(size(x,1)), tmpadd);
	addtable = [tmpaddmeans; tmpaddstderrs];
	addtable.Properties.RowNames = {'mean','stderr'};
	addtable.Properties.VariableNames = tmpadd.Properties.VariableNames;
	
	[~,adddispnames,~] = simul.get_simseries_array(addtable);
	[~,keep] = intersect( addtable.Properties.VariableNames, adddispnames, 'stable');
	N_keep = length(keep);
	
	statsout_exog{j}{end+(1:N_keep),:}=nan;
	statsout_exog{j}{end-N_keep+1:end,1}=addtable{'mean',keep}';
	statsout_exog{j}.Properties.RowNames(end-N_keep+1:end) = adddispnames;
	
	stderrs_exog{j}{end+(1:N_keep),:}=nan;
	stderrs_exog{j}{end-N_keep+1:end,1}=addtable{'stderr',keep}';
	stderrs_exog{j}.Properties.RowNames(end-N_keep+1:end) = adddispnames;
	
	addtables{j} = addtable;
end

freq_means = mean(freq_exog,2);