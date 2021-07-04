function print_iteration(iter, local_err, rs, f_evals, debug_errs)
%Print the information of each ftt iteration
%
%Tiangang Cui, August, 2019

fprintf('als=%2d, max_local_error=%3.3e, mean_local_error=%3.3e, max_rank=%d, cum#fevals=%3.3e\n', ...
    iter, max(local_err), mean(local_err), max(rs), f_evals);
if ~isempty(debug_errs)
    fprintf('als=%2d, max_debug_error=%3.3e, mean_debug_error=%3.3e\n', ...
        iter, debug_errs(1), debug_errs(2));
end
fprintf('\n');

end