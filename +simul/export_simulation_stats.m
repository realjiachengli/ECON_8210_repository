function out = export_simulation_stats(params,tabout,errtab,output_dir,resfile,prefix)

if nargin<6
	prefix='statsexog';
end

if ~strcmp( output_dir(end), '/') && ~strcmp( output_dir(end), '\')
	output_dir = [output_dir, '/'];
end

outstats=[output_dir,prefix,'_',resfile];
errstats=[output_dir,'errstats_',resfile];

for j=1:numel(tabout)
%      csvoutname=[outstats_exog,'_',num2str(j),'.csv'];
	writetable(tabout{j},outstats,'WriteRowNames',1,'FileType','spreadsheet','Sheet',j);
end
%writetable(paramout,[outstats_exog,'_params.csv'],'WriteRowNames',1);
writetable(params,outstats,'WriteRowNames',1,'FileType','spreadsheet','Sheet','params');
%writetable(errtab,errstats);
if ~isempty(errtab)
	writetable(errtab,errstats,'FileType','spreadsheet');
end