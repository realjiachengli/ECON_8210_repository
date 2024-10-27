function [slice,colnames] = estimate_sample_moments(simtmp,idx)

nvars = size(simtmp,2);
N_stats = 19;
colnames={'mean','std','skew','corrG','corrG_1','corrY','corrY_1','AC', 'min_obs', 'p01_obs', 'p05_obs', 'p10_obs', 'p25_obs', 'p50_obs', 'p75_obs', 'p90_obs', 'p95_obs', 'p99_obs', 'max_obs'};

slice = nan(nvars,N_stats);

if isempty(simtmp)
	return;
end
	
% fill with stats
slice(:,1)=mean(simtmp,1,'omitnan')';
slice(:,2)=std(simtmp,[],1,'omitnan')';
slice(:,3)=skewness(simtmp,[],1)';
% contemp and first-order autocorrelations
if size(simtmp,1)>2
	autocorrm=corrcoef([simtmp(2:end,:),simtmp(1:end-1,:)]);
	conm=autocorrm(1:nvars,1:nvars);
	lagm=autocorrm(nvars+1:end,1:nvars);
else
	conm=nan(nvars,nvars);
	lagm=nan(nvars,nvars);
end
% corr with shocks
slice(:,4:5)=[conm(:,idx(1)),lagm(idx(1),:)'];
% corr with Y
slice(:,6:7)=[conm(:,idx(2)),lagm(idx(2),:)'];
% vector with fo autocorr
slice(:,8)=diag(lagm);
%percentiles
if size(simtmp,1)>1
	slice(:,9:19)=prctile(simtmp,[0,1,5,10,25,50,75,90,95,99,100])';
else
	slice(:,9:19)=repmat(simtmp',1,11);
end

end