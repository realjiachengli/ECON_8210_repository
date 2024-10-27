k_next = [1,2]';
inv_next = [3,4]';
Vf = @(a) a(:, 1);


num_exog_states = 3;

% Create a vector for exogenous states (1 to 9)
exog_states = 1:num_exog_states;

% Use ndgrid to generate all combinations of exog_states, k_next, and inv_next
[ExogGrid, KGrid, InvGrid] = ndgrid(exog_states, k_next, inv_next);

% Reshape these grids into an (npts*num_exog_states) by 3 matrix for evaluation
% We need to transpose grids since ndgrid arranges the first dimension (exog_states)
future_states = [ExogGrid(:), KGrid(:), InvGrid(:)];

% Evaluate the value function Vf at all future states
future_values_vector = Vf(future_states);

% Reshape the result back into an npts by 9 matrix (one column per exogenous state)
future_values = reshape(future_values_vector, 2, num_exog_states);