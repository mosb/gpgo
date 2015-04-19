inproblems = [];
%branin
problem.name = 'branin';
inproblems = [inproblems; problem];

%camelback
problem.name = 'camelback';
inproblems = [inproblems; problem];

%goldstein price
problem.name = 'goldsteinprice';
inproblems = [inproblems; problem];

%har 3
problem.name = 'hartman3';
inproblems = [inproblems; problem];

%har 6
problem.name = 'hartman6';
inproblems = [inproblems; problem];

%shekel 10
problem.name = 'shekel10';
inproblems = [inproblems; problem];

%shekel 7
problem.name = 'shekel7';
inproblems = [inproblems; problem];

%shekel 5
problem.name = 'shekel5';
inproblems = [inproblems; problem];

%shubert
problem.name = 'shubert';
inproblems = [inproblems; problem];

%griewank
problem.name = 'griewank2';
inproblems = [inproblems; problem];

%griewank
problem.name = 'griewank5';
inproblems = [inproblems; problem];

%griewank
problem.name = 'griewank10';
inproblems = [inproblems; problem];

%ackley 
problem.name = 'ackley2';
inproblems = [inproblems; problem];

%ackley 
problem.name = 'ackley5';
inproblems = [inproblems; problem];

%ackley 
problem.name = 'ackley10';
inproblems = [inproblems; problem];

%rast 2
problem.name = 'rast';
inproblems = [inproblems; problem];

inproblems = inproblems(probinds);

fid = fopen(['usresults-1-1-25-' num2str(min(probinds)) ...
             '-' num2str(max(probinds))],'w');  
for maxevals = [10 20 50]
  disp(['maxevals: ' num2str(maxevals)]);
  fprintf(fid,['maxevals: ' num2str(maxevals) '\n']);
  for pind=1:numel(inproblems)
    inproblem = inproblems(pind);
    eval(['load tomlabresults/out-' inproblem.name '-' num2str(maxevals)]);
    thisproblems = outproblems;
    disp(['problem: ' inproblem.name]);
    outproblems = [];
    for off = 1:25
      display(['offset number: ' num2str(off)]);
      problem = thisproblems(off);

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
      
      problem.oursamples = [1];
      covvys = [struct('covfn',@(x)ndimsqdexp_isotropic_cov_fn(x,problem.n)); ...
        struct('covfn',@(x)ndimsqdexpperiodic_isotropic_cov_fn(x,problem.n))];

      for samples = problem.oursamples;
        for covvyind = 1:length(covvys)
          disp(['samples: ' num2str(samples)]);
          problem.covvy = covvys(covvyind);
          problem.covvy

          options = struct( ...
              'MaxFunEvals', maxevals * problem.n, ...
              'ValueTol', 1e-3, ...
              'UseNAG', 0, ...
              'Steps', 1 ...
              );
            
%           if (strcmp(problem.name, 'shubert'))
%             problem.outscale = 50;
%           end
% 
%           if (strfind(problem.name, 'griewank') > 0)
%             problem.outscale = 30000;
%           end
  
          problem.outscale = 1e8;

          for inscale = [5 10 20]
            disp(['inscale: ' num2str(inscale)]);
            problem.inscale = max(problem.upprbnd - problem.lowrbnd) / inscale;
            %problem.instd = max(log(problem.inscale), 1);
            problem.instd = 1;
            problem.outstd = 3;
            
            problem.covvy.hyperparams(1) = ...
              struct('name','logInputScale','priorMean',log(problem.inscale),'priorSD',problem.instd,'NSamples',samples,'type','real');
            problem.covvy.hyperparams(2) = ...
              struct('name','logOutputScale','priorMean',log(problem.outscale),'priorSD',problem.outstd,'NSamples',samples,'type','real');
            problem.covvy.hyperparams(3) = ...
              struct('name','mean','priorMean',0,'priorSD',3,'NSamples',1,'type','real');
            problem.covvy.hyperparams(4) = ...
              struct('name','logNoiseSD','priorMean',log(0.0001),'priorSD',0.1,'NSamples',1,'type','real');

            tic
            fn=problem.f;
            x0=problem.x0;
            lowrbnd = problem.lowrbnd;
            upprbnd = problem.upprbnd;
            covvy=problem.covvy;
            plotit = 0;
            maximise;
            o = outputs;
            
            problem.ouroutputs = [problem.ouroutputs; o];
            problem.ourresults = [problem.ourresults; -max(o.ys)]; 
            problem.ourerrors = [problem.ourerrors; ...
              abs((problem.f(problem.x0) - problem.ourresults(end)) / (problem.f(problem.x0) - problem.optimum))];
            disp(['our error: ' num2str(problem.ourerrors(end))]);
            problem.ourmxresults = [problem.ourmxresults; fn(o.mxLoc)]; 
            problem.ourmxerrors = [problem.ourmxerrors; ...
              abs((problem.f(problem.x0) - problem.ourmxresults(end)) / (problem.f(problem.x0) - problem.optimum))];
            disp(['our mxerror: ' num2str(problem.ourmxerrors(end))]);
            toc
          end
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

      outstring = [problem.name ' usmx:'];
      for nind = 1:numel(problem.ourmxerrors)
        outstring = [outstring ' ' num2str(problem.ourmxerrors(nind))];
      end                                
      disp(outstring);
      fprintf(fid, [outstring '\n']);

      outstring = [problem.name ' us mx best: ' num2str(max(problem.ourmxerrors))];
      disp(outstring);
      fprintf(fid, [outstring '\n']);
      
      outproblems = [outproblems; problem];
    end
    cd ourresults
    eval(['save out-1-1-25-' problem.name '-' num2str(maxevals) ' outproblems']);
    cd ..
  end
end
