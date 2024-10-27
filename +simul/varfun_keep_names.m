function functab = varfun_keep_names(func,tab)
	functab = varfun(func,tab,'GroupingVariables','per');
	functab.GroupCount=[];
		
	renamevars = setdiff(tab.Properties.VariableNames, 'per','stable');
	functab.Properties.VariableNames(2:end) = renamevars;
	functab.run=[];
end
