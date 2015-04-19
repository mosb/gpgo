inproblems = [];

%branin
problem.name = 'branin';
problem.periodic = false;
inproblems = [inproblems; problem];
%camelback
problem.name = 'camelback';
problem.periodic = false;
inproblems = [inproblems; problem];
%goldstein price
problem.name = 'goldsteinprice';
problem.periodic = false;
inproblems = [inproblems; problem];
%har 3
problem.name = 'hartman3';
problem.periodic = false;
inproblems = [inproblems; problem];
%har 6
problem.name = 'hartman6';
problem.periodic = false;
inproblems = [inproblems; problem];
%shekel 10
problem.name = 'shekel10';
problem.periodic = false;
inproblems = [inproblems; problem];
%shekel 7
problem.name = 'shekel7';
problem.periodic = false;
inproblems = [inproblems; problem];
%shekel 5
problem.name = 'shekel5';
problem.periodic = false;
inproblems = [inproblems; problem];
%shubert
problem.name = 'shubert';
problem.periodic = false;
inproblems = [inproblems; problem];
%griewank
problem.name = 'griewank2';
problem.periodic = false;
inproblems = [inproblems; problem];
%griewank
problem.name = 'griewank5';
problem.periodic = false;
inproblems = [inproblems; problem];
%griewank
%problem.name = 'griewank10';
%inproblems = [inproblems; problem];
%ackley 
problem.name = 'ackley2';
problem.periodic = false;
inproblems = [inproblems; problem];
%ackley 
problem.name = 'ackley5';
problem.periodic = false;
inproblems = [inproblems; problem];
%ackley 
%problem.name = 'ackley10';
%inproblems = [inproblems; problem];
%rast 2
problem.name = 'rast';
problem.periodic = false;
inproblems = [inproblems; problem];

inproblems = inproblems(probinds);

for maxevals = [10]
  fid = fopen(['mult-usresults-1-10-' num2str(probinds(1)) ...
               '-' num2str(probinds(end)) '-' num2str(maxevals)],'w');  
  disp(['maxevals: ' num2str(maxevals)]);
  todo = 1;
  while (todo > 0)
    todo = 0;
    for pind=1:numel(inproblems)
      inproblem = inproblems(pind);
      eval(['load tomlabresults/out-' inproblem.name '-' num2str(maxevals)]);
      thisproblems = outproblems;
			for i=1:numel(thisproblems)
				thisproblems(i).periodic = inproblem.periodic;
			end
      disp(['problem: ' inproblem.name]);

      cd ourresults
      if (exist(['mult-out-1-10-' inproblem.name '-' num2str(maxevals) '.mat']) == 2)
        eval(['load mult-out-1-10-' inproblem.name '-' num2str(maxevals)]);
      else
        outproblems = [];
        off = 1;
      end
      cd ..
      
      if (off > 10)
        continue
      end
      
      display(['offset number: ' num2str(off)]);
      problem = thisproblems(off);

     	problem.forig = problem.f;
    	problem.f = @(x)-problem.forig(x);

      problem.egoerror = ...
        (problem.forig(problem.x0) - problem.egoresults.f_k) / ...
					(problem.forig(problem.x0) - problem.optimum);
      problem.rbferror = ...
        (problem.forig(problem.x0) - problem.rbfresults.f_k) / ...
					(problem.forig(problem.x0) -problem.optimum);
      problem.directerror = ...
        (problem.forig(problem.x0) - problem.directresults.f_k) / ...
					(problem.forig(problem.x0) - problem.optimum);
			
			problem.egoerror = max(min(problem.egoerror,1),0);
			problem.rbferror = max(min(problem.rbferror,1),0);
			problem.directerror = max(min(problem.directerror,1),0);

      outstring = [problem.name ' ego: ' num2str(problem.egoerror) ...
                                ' rbf: ' num2str(problem.rbferror) ... 
                                ' direct: ' num2str(problem.directerror)];
      disp(outstring);
      fprintf(fid, [outstring '\n']);

      if (~isfield(problem,'ourerrors')); problem.ourerrors = []; end
      if (~isfield(problem,'ourmxerrors')); problem.ourmxerrors = []; end
      if (~isfield(problem,'ourresults')); problem.ourresults = []; end
      if (~isfield(problem,'ourrmxesults')); problem.ourmxresults = []; end
      if (~isfield(problem,'ouroutputs')); problem.ouroutputs = []; end

      problem.ourinsamples = [3];
      problem.ouroutsamples = [5];
  		covvys = [struct('covfn',@(x)ndimsqdexp_isotropic_cov_fn(x,problem.n))];
			
      for sind = 1:numel(problem.ourinsamples)
        for covvyind = 1:length(covvys)
          derivobs = false;
          insamples = problem.ourinsamples(sind);
          outsamples = problem.ouroutsamples(sind);
          disp(['in samples: ' num2str(insamples) ...
								' out samples: ' num2str(outsamples) ...
								' deriv obs: ' num2str(derivobs)]);

          problem.covvy = covvys(covvyind);
          problem.covvy

          options = struct( ...
              'MaxFunEvals', maxevals * problem.n, ...
              'ValueTol', 1e-3, ...
              'UseNAG', 0, ...
              'Steps', 2, ...
              'ForEach', 3, ...
              'DerivObs', derivobs ...
              );

          problem.outscale = exp(6);
          problem.outstd = 4;

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
            struct('name','mean', ...
                   'priorMean',0, ...
                   'priorSD',3, ...
                   'NSamples',1, ...
                   'type','real');
          problem.covvy.hyperparams(4) = ...
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
          maximise;
          o = outputs;
          %o =
          %maximise(problem.f,problem.x0,problem.lowrbnd,problem.upprbnd,problem.covvy,options);
					if (o.successful == false)
						disp('failed!');
						fprintf(fid,'failed!\n');
					end

          problem.ouroutputs = [problem.ouroutputs; o];
          problem.ourresults = [problem.ourresults; -max(o.ys)]; 
          problem.ourerrors = [problem.ourerrors; ...
            (problem.forig(problem.x0) - problem.ourresults(end)) / ...
										(problem.forig(problem.x0) - problem.optimum)];
					problem.ourerrors = max(min(problem.ourerrors,1),0);
          disp(['our error: ' num2str(problem.ourerrors(end))]);
          problem.ourmxresults = [problem.ourmxresults; -problem.f(o.mxLoc)]; 
          problem.ourmxerrors = [problem.ourmxerrors; ...
            (problem.forig(problem.x0) - problem.ourmxresults(end)) / ...
										(problem.forig(problem.x0) - problem.optimum)];
					problem.ourmxerrors = max(min(problem.ourmxerrors,1),0);
          disp(['our mxerror: ' num2str(problem.ourmxerrors(end))]);
          toc
        end
      end
      outstring = [problem.name ' ego: ' num2str(problem.egoerror) ...
									 ' rbf: ' num2str(problem.rbferror) ... 
									 ' direct: ' num2str(problem.directerror)];
      disp(outstring);
			
      outstring = [problem.name ' us:'];
      for nind = 1:numel(problem.ourerrors)
        outstring = [outstring ' ' num2str(problem.ourerrors(nind))];
      end                                
      disp(outstring);
      fprintf(fid, [outstring '\n']);
			
      outstring = [problem.name ' us best: ' num2str(max(problem.ourerrors))];
      disp(outstring);
      fprintf(fid, [outstring '\n']);
			
      outproblems = [outproblems; problem];
      off = off + 1;
      if (off < 11)
        todo = 1;
      end
      cd ourresults
      eval(['save -v7.3 mult-out-1-10-' problem.name '-' num2str(maxevals) ' outproblems off']);
      cd '..'
		end
					
	end
	fclose(fid);
end
