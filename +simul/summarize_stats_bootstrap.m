function [statsout_exog, stderrs_exog, dispnames, freq_means, addtables] = summarize_stats_bootstrap( statscell, freqcell, addcell )

N_runs = numel(statscell);
N_subsamples = numel(statscell{1});

blanktable = statscell{1}{1};
sz = size(blanktable);
N_vars = sz(1);
N_stats = sz(2);

statsout_exog = cell(N_subsamples,1);
stderrs_exog = cell(N_subsamples,1);
addtables = cell(N_subsamples,1);



for j=1:N_subsamples 
	statstmp=zeros(N_vars,N_stats,N_runs);
	tmp_addcell = cellfun(@(c)c{j}, addcell, 'UniformOutput', false);
	%for ii=1:N_runs
	parfor ii=1:N_runs
		statstmp(:,:,ii) = statscell{ii}{j}{:,:};
	end
	
	tmpmeans=mean(statstmp,3,'omitnan');
	statsout_exog{j} = blanktable;
	statsout_exog{j}{:,:} = tmpmeans;
	
	tmpstderrs=std(statstmp,[],3,'omitnan') / sqrt(N_runs);
	stderrs_exog{j} = blanktable;
	stderrs_exog{j}{:,:} = tmpstderrs;

	% Deal with additional stats
	tmpadd = vertcat( tmp_addcell{:} );
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

freq_means = mean([freqcell{:}],2);

dispnames = statsout_exog{1}.Properties.RowNames;