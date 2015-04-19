inproblems = [];

%shubert
problem.name = 'shubert';
problem.periodic = true;
inproblems = [inproblems; problem];
%griewank
problem.name = 'griewank2';
problem.periodic = true;
inproblems = [inproblems; problem];
%ackley 
problem.name = 'ackley2';
problem.periodic = true;
inproblems = [inproblems; problem];
%rast 2
problem.name = 'rast';
problem.periodic = true;
inproblems = [inproblems; problem];
%ackley 
problem.name = 'ackley5';
problem.periodic = true;
inproblems = [inproblems; problem];
%griewank
problem.name = 'griewank5';
problem.periodic = true;
inproblems = [inproblems; problem];

inproblems = inproblems(probinds);

for maxevals = [10]
  disp(['maxevals: ' num2str(maxevals)]);
	for pind=1:numel(inproblems)
		inproblem = inproblems(pind);
		eval(['load ourresults/fixed-out-1-10-' inproblem.name '-' num2str(maxevals)]);
		thisproblems = outproblems;
		disp(['problem: ' inproblem.name]);
		
		for offind = 8:10
      display(['offset number: ' num2str(offind)]);
      problem = thisproblems(offind);
			
      outstring = [problem.name ' ego: ' num2str(problem.egoerror) ...
									 ' rbf: ' num2str(problem.rbferror) ... 
									 ' direct: ' num2str(problem.directerror)];
      disp(outstring);
			
      problem.ourinsamples = [3];
      problem.ouroutsamples = [3];
			covvys = [struct('covfn',@(x)ndimsqdexp_isotropic_cov_fn(x,problem.n)); ...
								struct('covfn',@(x)ndimsqdexpperiodic_isotropic_cov_fn(x,problem.n))]; 
			
      for sind = 1:numel(problem.ourinsamples)
        for covvyind = 2
          insamples = problem.ourinsamples(sind);
          outsamples = problem.ouroutsamples(sind);
          disp(['in samples: ' num2str(insamples) ...
								' out samples: ' num2str(outsamples)]);
					
          problem.covvy = covvys(covvyind);
          problem.covvy
					
          options = struct( ...
              'MaxFunEvals', maxevals * problem.n, ...
              'ValueTol', 1e-3, ...
              'UseNAG', 0, ...
              'Steps', 1, ...
							'ForEach', 10, ...
              'DerivObs', false ...
              );
					
			    problem.outscale = exp(2);
          problem.outstd = 1;
					
          problem.inscale = max(problem.upprbnd - problem.lowrbnd) / 10;
          problem.instd = 0.75;
					
					problem.covvy.hyperparams(1) = ...
							struct('name','logInputScale', ...
										 'priorMean',log(problem.inscale), ...
										 'priorSD',problem.instd, ...
										 'NSamples',insamples, ...
										 'type','real');
					problem.covvy.hyperparams(2) = ...
							struct('name','logOutputScale', ...
										 'priorMean',log(problem.outscale), ...
										 'priorSD',problem.outstd, ...
										 'NSamples',outsamples, ...
										 'type','real');
					problem.covvy.hyperparams(3) = ...
							struct('name','logInputScale2', ...
										 'priorMean',log(problem.inscale), ...
										 'priorSD',problem.instd, ...
										 'NSamples',insamples, ...
										 'type','real');
					problem.covvy.hyperparams(4) = ...
							struct('name','logOutputScale2', ...
										 'priorMean',log(problem.outscale), ...
										 'priorSD',problem.outstd, ...
										 'NSamples',outsamples, ...
										 'type','real');
					problem.covvy.hyperparams(5) = ...
							struct('name','mean', ...
										 'priorMean',0, ...
										 'priorSD',3, ...
										 'NSamples',1, ...
										 'type','real');
					problem.covvy.hyperparams(6) = ...
							struct('name','logNoiseSD', ...
										 'priorMean',log(0.00001), ...
										 'priorSD',0.1, ...
										 'NSamples',1, ...
										 'type','real');
					
          tic
          fn=problem.f;
          x0=problem.x0;
          lowrbnd = problem.lowrbnd;
          upprbnd = problem.upprbnd;
          covvy=problem.covvy;
          plotit = 0;
          maximise_direct;
          o = outputs;
					if (o.successful == false)
						disp('failed!');
					end
					
					ind = 4 + sind;
					
          problem.ouroutputs(ind) = o;
          problem.ourresults(ind) = -max(o.ys); 
          problem.ourerrors(ind) = ...
										(problem.forig(problem.x0) - problem.ourresults(ind)) / ...
										(problem.forig(problem.x0) - problem.optimum);
					problem.ourerrors = max(min(problem.ourerrors,1),0);
          disp(['our error: ' num2str(problem.ourerrors(ind))]);
          problem.ourmxresults(ind) = -problem.f(o.mxLoc); 
          problem.ourmxerrors(ind) = ...
										(problem.forig(problem.x0) - problem.ourmxresults(ind)) / ...
										(problem.forig(problem.x0) - problem.optimum);
					problem.ourmxerrors = max(min(problem.ourmxerrors,1),0);
          disp(['our mxerror: ' num2str(problem.ourmxerrors(ind))]);
          toc
        end
      end

      outstring = [problem.name ' us:'];
      for nind = 1:numel(problem.ourerrors)
        outstring = [outstring ' ' num2str(problem.ourerrors(nind))];
      end                                
      disp(outstring);
			
      outstring = [problem.name ' us best: ' num2str(max(problem.ourerrors))];
      disp(outstring);
			
      outstring = [problem.name ' usmx:'];
      for nind = 1:numel(problem.ourmxerrors)
        outstring = [outstring ' ' num2str(problem.ourmxerrors(nind))];
      end                                
      disp(outstring);
			
      outstring = [problem.name ' us mx best: ' num2str(max(problem.ourmxerrors))];
      disp(outstring);
			
      outproblems(offind)= problem;
      cd ourresults
      eval(['save periodicfix-out-1-10-' problem.name '-' num2str(maxevals) ' outproblems off']);
      cd '..'
		end
	end
end
