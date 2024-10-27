function [simarray,dispnames,nvars] = get_simseries_array(simseries)

single_col = varfun(@(x)size(x,2)==1,simseries,'OutputFormat','uniform');
dispnames=simseries.Properties.VariableNames(single_col);
nvars = length(dispnames);
simarray = simseries{:,single_col};