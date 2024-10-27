function paths = export_paths(tab,varnames,NT_sim,N_runs)
	paths=struct;
    for v=1:length(varnames)
        thisvar=zeros(N_runs,NT_sim);
        for p=1:N_runs
            thisvar(p,:)=tab{tab.run==p,varnames{v}};
        end
        paths.(varnames{v})=thisvar;
    end
end
