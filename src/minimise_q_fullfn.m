function qs = minimise_q_fullfn(covvy,samples,x,x_data,y_data,evaluation,mx)

num_hypersamples = numel(covvy.hypersamples);
qs = nan(num_hypersamples,1);

if isempty(samples)
    return
end

for sample = samples
  cholK = downdatechol(covvy.hypersamples(sample).cholK, evaluation);
  covvy.hypersamples(sample).cholK = cholK;
  covvy.hypersamples(sample).datatwothirds = solve_chol(cholK, y_data(1:evaluation - 1));
  
  qs(sample) = max(eps, -negvalue(x, x_data(1:evaluation - 1, :), y_data(1:evaluation - 1), covvy, sample, mx) - mx);
        % I want qs to be expected improvement so that it is strictly
        % positive, ROMAN can you check that I have done it right.
end