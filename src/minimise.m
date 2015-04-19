%minimise function fn by sequentially greedily selecting the most valuable
%observation according to a GP fit to fn. 

function outputs = minimise(fn, x_0, lower_bound, upper_bound, covvy, options)
  
  fn = @(x)-fn(x);
  dim = length(x_0);
  
  default_options = struct('function_evaluations', 10 * dim, ...
                           'lookahead_steps', 1, ...
                           'derivative_observations', true, ...
                           'input_scales', 1, 'plot', true );
 
  names = fieldnames(default_options);
  for i = 1:length(names);
    name = names{i};
    if (~isfield(options, name))
      options.(name) = default_options.(name);
    end
  end
  
  x = x_0;
  
  if (options.derivative_observations)
    derivative_directions = zeros(options.function_evaluations, dim);
  end
  
  augment = @(x) [x zeros(options.derivative_observations * size(x))];
  
  x_locations = zeros(options.function_evaluations, dim);
  y_data = zeros(options.function_evaluations, 1);

  actual_x_locations = zeros(options.function_evaluations, dim);
  actual_y_data = zeros(options.function_evaluations, 1);
  
  distances = [];

  % actual_y_data differs from y_data in that where we replace an
  % observation with a derivative observation, the original
  % observation is stored in actual_y_data, where the derivative obs
  % goes in y_data. The location of the actual observation goes in
  % actual_x_locations.
  
  actual_observation_flags = false(options.function_evaluations, 1);

  covvy = hyperparams(covvy);

  length_scales = ones(1, dim);
  input_scale_flags = cellfun(@(x) strcmp(x, 'logInputScale'), covvy.names);
  if any(input_scale_flags)
    length_scales(:) = mean(covvy.hyperparams(find(input_scale_flags)).samples);
  else
    for i = 1:dim
      this_input_scale = find(cellfun(@(x) strcmp(x, ['logInputScale' ...
                          num2str(i)]), covvy.names));
      length_scales(i) = mean(covvy.hyperparams(this_input_scale).samples);
    end
  end
  personal_space = options.input_scales * exp(length_scales);
  
  covvy = bmcparams(covvy);
  
  for evaluation = 1:options.function_evaluations
%    tic
%    disp(['evaluation: ' num2str(evaluation)]);
    y = fn(x);
    
    % find the matrix of distances of every observation from every
    % other regular observation
    newdists = ((repmat(x, evaluation - 1, 1) - ...
                 actual_x_locations(1:evaluation - 1, :)).^2 * ...
                personal_space'.^-2);
    
    actual_x_locations(evaluation, :) = x;
    actual_y_data(evaluation) = y;
    
    if (options.derivative_observations)
      % Is this new observation is too close to any existing function
      % observation? Note this means we might still get problems due
      % to derivative observations getting too close, but these are
      % much better behaved.
      if any(newdists(actual_observation_flags) < 1)
        % We are going to replace this observation with a derivative
        % observation
        [minimum, index] = min(newdists);
        direction = (x - actual_x_locations(index, :));
        x = actual_x_locations(index, :) + 0.5 * direction; % at the midpoint
        y = (y - actual_y_data(index)) / sqrt(direction.^2 * ones(dim, 1));
        actual_observation_flags(evaluation) = false;
      else
        direction = zeros(1, dim);
        actual_observation_flags(evaluation) = true;
      end
      derivative_directions(evaluation, :) = direction;
    end
    
    distances = [distances, newdists; newdists', nan];

    x_locations(evaluation, :) = x;
    y_data(evaluation) = y;

    if (options.derivative_observations)
      x_data = [x_locations, derivative_directions];
    else
      x_data = x_locations;
    end
    
    covvy = gpparams(x_data(1:evaluation, :), y_data(1:evaluation), ...
                     covvy, 'update', evaluation);

    rho = weights(covvy);

    [maximum, index] = max(actual_y_data(1:evaluation));
    maximum_location = actual_x_locations(index, :);
    
    if (evaluation == options.function_evaluations)
      continue
    elseif (evaluation == 1) % we should just go to a corner
			x = lower_bound;
			continue
		end
	
		if (options.lookahead_steps == 2)
      Prob.user.lowrbnd = lower_bound;
      Prob.user.upprbnd = upper_bound;
      Prob.user.num_samples = 7;
    end

    Prob.user.rho = rho;
    Prob.user.XData = x_data(1:evaluation, :);
    Prob.user.YData = y_data(1:evaluation);
    Prob.user.covvy = covvy;
    Prob.user.mx = maximum;

  	if (options.lookahead_steps == 1)
			Problem.f = @(XStar)weightednegval_direct(augment(XStar)', Prob);
		else
			Problem.f = @(XStar)multi_step_negvalue_direct_2(augment(XStar)', Prob);
    end

    opts.maxevals = 10000;
		opts.showits = 0;
		bounds = [lower_bound; upper_bound]';
		
    [retmin, minimum_location] = Direct(Problem, bounds, opts);
		x = minimum_location';
    
    if (dim ==1 && options.plot)
			minimise_val_plot;
    elseif (dim == 2 && options.plot && mod(evaluation, 10) == 0)
			plotstuff;
    end
  
%  toc
  end

outputs.xs = x_data;
outputs.ys = -y_data;
outputs.covvy = covvy;
outputs.min = -maximum;
outputs.minLoc = maximum_location;
outputs.evaluation = evaluation;
outputs.successful = true;
