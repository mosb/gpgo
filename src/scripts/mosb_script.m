%clear
%close all
%drawnow
%cd '/raidy/work/common/tomlab'
%try
%  startup
%catch
%end
%cd '/raidy/work/optimization'

probinds=1

for maxevals = [25]
  disp(['maxevals: ' num2str(maxevals)]);
  inproblems = [];
  outproblems = [];

  % bra.m           Branin
  % cam.m           Four-hump camel
  % gpr.m           Goldstein-Price
  % hm3.m           Hartman3
  % hm6.m           Hartman6
  % ros.m           Rosenbrock
  % s10.m           Shekel10
  % sh5.m           Shekel5
  % sh7.m           Shekel7
  % shu.m           Shubert
  
  % 1: branin
  % 3.9416
  % 1.9988
  % 2: camelback
  % 10.3817
  % 9.9083
  % 3: goldsteinprice
  % 16.2883
  % 16.0795
  % 4: hartman3
  % -0.028464
  % -2.6555
  % 5: hartman6
  % -1.261
  % -3.0539
  % 6: shekel10
  % -1.8267
  % -3.4975
  % 7: shekel7
  % -2.0855
  % -3.6105
  % 8: shekel5
  % -2.3692
  % -4.0223
  % 9: shubert
  % 3.2952
  % 0.74088
  % 10: griewank2
  % 12.2031
  % 10.9883
  % 11: griewank5
  % 11.7449
  % 9.7765
  % 12: griewank10
  % 11.375
  % 9.0873
  % 13: ackley2
  % 0.83342
  % -1.3736
  % 14: ackley5
  % -0.35828
  % -2.2739
  % 15: ackley10
  % -1.1219
  % -3.1675
  % 16: rast
  % 3.4336
  % 1.8496

  %michalewics
  %problem.shortname = 'michalewics 5';
  %problem.n = 5;
  %problem.lowrbnd = zeros(1,problem.n);
  %problem.upprbnd = pi * ones(1,problem.n);
  %problem.f = @mich;
  %problem.optimum = -4.687658;
  %problems = [problems; problem];

  %branin
  problem.name = 'branin';
  problem.outscale = 5;
  problem.outstd = 2;
  inproblems = [inproblems; problem];
  
  %camelback
  problem.name = 'camelback';
  problem.outscale = 10;
  problem.outstd = 2;
  inproblems = [inproblems; problem];
  
  %goldstein price
  problem.name = 'goldsteinprice';
  problem.outscale = 10;
  problem.outstd = 2;
  inproblems = [inproblems; problem];
  
  %har 3
  problem.name = 'hartman3';
  problem.outscale = 1;
  problem.outstd = 1;
  inproblems = [inproblems; problem];

  %har 6
  problem.name = 'hartman6';
  problem.outscale = 1;
  problem.outstd = 1;
  inproblems = [inproblems; problem];
 
  %shekel 10
  problem.name = 'shekel10';
  problem.outscale = 1;
  problem.outstd = 1;
  inproblems = [inproblems; problem];

  %shekel 7
  problem.name = 'shekel7';
  problem.outscale = 1;
  problem.outstd = 1;
  inproblems = [inproblems; problem];
  
  %shekel 5
  problem.name = 'shekel5';
  problem.outscale = 1;
  problem.outstd = 1;
  inproblems = [inproblems; problem];

  %shubert
  problem.name = 'shubert';
  problem.outscale = 5;
  problem.outstd = 2;
  inproblems = [inproblems; problem];

  %griewank
  problem.name = 'griewank2';
  problem.outscale = 10;
  problem.outstd = 2;
  inproblems = [inproblems; problem];

  %griewank
  problem.name = 'griewank5';
  problem.outscale = 10;
  problem.outstd = 2;
  inproblems = [inproblems; problem];

  %griewank
  problem.name = 'griewank10';
  problem.outscale = 10;
  problem.outstd = 2;
  inproblems = [inproblems; problem];
  
  %ackley 
  problem.name = 'ackley2';
  problem.outscale = 1;
  problem.outstd = 1;
  inproblems = [inproblems; problem];

  %ackley 
  problem.name = 'ackley5';
  problem.outscale = 1;
  problem.outstd = 1;
  inproblems = [inproblems; problem];

  %ackley 
  problem.name = 'ackley10';
  problem.outscale = 1;
  problem.outstd = 1;
  inproblems = [inproblems; problem];
  
  %rast 2
  problem.name = 'rast';
  problem.outscale = 5;
  problem.outstd = 2;
  inproblems = [inproblems; problem];
  
  inproblems = inproblems(probinds);

%   for pind = 1:numel(problems)
%     thisproblem = [];
%     for off = 1:50
%       display(['offset number: ' num2str(off)]);
%       problem = problems(pind);
%       disp(['problem: ' problem.name]);
% 
%       problem.offset = zeros(1,problem.n);
%       for i = 1:problem.n
%         range = problem.canAdd(i) + problem.canSubtract(i);
%         num = rand * range;
%         if (num > problem.canAdd(i))
%           num = -(num - problem.canAdd(i));
%         end
%         problem.offset(i) = num;
%       end
% 
%       problem.lowrbnd = problem.lowrbnd + problem.offset;
%       problem.upprbnd = problem.upprbnd + problem.offset;
%       problem.bounds = [problem.lowrbnd; problem.upprbnd]';
%       problem.x0 = (problem.upprbnd - problem.lowrbnd) / 2 + problem.lowrbnd;
%       thisproblem = [thisproblem; problem];
%     end
%     eval(['save ' problem.name ' thisproblem']);
%   end
    
  fid = fopen('results','w');  

  for off = 1:50
    display(['offset number: ' num2str(off)]);
    for pind=1:numel(inproblems)
      inproblem = inproblems(pind);
      eval(['load testfunctions/' inproblem.name]);
      problem = thisproblem(off);
      problem.maxevals = maxevals;
      disp(['problem: ' problem.name]);

      problem.inscale = max(problem.upprbnd - problem.lowrbnd) / 10;
      problem.instd = max(log(problem.inscale) * 2, 1);
      problem.outscale = inproblem.outscale;
      problem.outstd = inproblem.outstd;
      
      problem.ourerrors = [];
      problem.ouroutputs = [];
      problem.oursamples = [3];
      covvys = [struct('covfn',@(x,varargin)ndimsqdexp_isotropic_cov_fn(x,problem.n,varargin{:})); ...
        struct('covfn',@(x,varargin)ndimsqdexpperiodic_isotropic_cov_fn(x,problem.n,varargin{:}))];

      for samples = problem.oursamples;
        for covvyind = 1:length(covvys)
          disp(['samples: ' num2str(samples)]);
          problem.covvy = covvys(covvyind);
          problem.covvy

          problem.covvy.hyperparams(1) = ...
            struct('name','logInputScale','priorMean',log(problem.inscale),'priorSD',problem.instd,'NSamples',samples,'type','real');
          problem.covvy.hyperparams(2) = ...
            struct('name','logOutputScale','priorMean',problem.outscale,'priorSD',problem.outstd,'NSamples',samples,'type','real');
          problem.covvy.hyperparams(3) = ...
            struct('name','mean','priorMean',0,'priorSD',3,'NSamples',1,'type','real');
          problem.covvy.hyperparams(4) = ...
            struct('name','logNoiseSD','priorMean',log(0.0001),'priorSD',0.1,'NSamples',1,'type','real');

          options = struct( ...
              'MaxFunEvals', maxevals * problem.n, ...
              'ValueTol', 1e-3, ...
              'UseNAG', 0, ...
              'Steps', 1 ...
              ); 
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
          problem.ourresult = -max(o.ys); 
          problem.ourerrors = [problem.ourerrors; ...
            abs((problem.f(problem.x0) - problem.ourresult) / (problem.f(problem.x0) - problem.optimum))];
          disp(['our error: ' num2str(problem.ourerrors(end))]);
          toc
        end
      end

%       Prob = glcAssign(problem.shortname,problem.lowrbnd,problem.upprbnd,problem.name,[],[],[],'',[],[],problem.x0);
%       Prob.optParam.MaxFunc = maxevals * problem.n;
%       Prob.optParam.IterPrint = 0;
%       
%       results = ego(Prob);
%       problem.egoresults = results;
%       problem.egoerror = ...
%           abs((problem.f(problem.x0) - results.f_k) / (problem.f(problem.x0) - problem.optimum));
%       disp(['ego error: ' num2str(problem.egoerror)]);
% 
%       results = rbfSolve(Prob);
%       problem.rbfresults = results;
%       problem.rbferror = ...
%           abs((problem.f(problem.x0) - results.f_k) / (problem.f(problem.x0) - problem.optimum));
%       disp(['ego error: ' num2str(problem.egoerror)]);
%  
%       results = glbSolve(Prob);
%       problem.directresults = results;
%       problem.directerror = ...
%           abs((problem.f(problem.x0) - results.f_k) / (problem.f(problem.x0) - problem.optimum));
%       disp(['direct error: ' num2str(problem.directerror)]);
% 
       outproblems = [outproblems; problem];
       for ind=1:numel(outproblems)
         outstring = [outproblems(ind).name ' us: '];
         for j=1:numel(outproblems(ind).ourerrors)
           outstring = [outstring num2str(outproblems(ind).ourerrors(j)) ' '];
         end
%         outstring = [outstring 'ego: ' num2str(outproblems(ind).egoerror) ...
%                                ' rbf: ' num2str(outproblems(ind).rbferror) ... 
%                                ' direct: ' num2str(outproblems(ind).directerror)];
         disp(outstring);
         fprintf(fid, [outstring '\n']);
       end
    end
  end
  %eval(['save results-' num2str(min(probinds)) '-' num2str(max(probinds)) ...
  %  '-' num2str(maxevals) ' outproblems']);
end

%       for i=1:problem.n
%         problem.covvy.hyperparams(i) = ...
%           struct('name',['logInputScale' num2str(i)],'priorMean',log(problem.inscale),'priorSD',problem.instd,'NSamples',samples,'type','real');
%       end
%       problem.covvy.hyperparams(problem.n+1) = ...
%         struct('name','logOutputScale','priorMean',log(problem.outscale),'priorSD',problem.outstd,'NSamples',samples,'type','real');
%       problem.covvy.hyperparams(problem.n+2) = ...
%         struct('name','mean','priorMean',0,'priorSD',3,'NSamples',1,'type','real');
%       problem.covvy.hyperparams(problem.n+3) = ...
%         struct('name','logNoiseSD','priorMean',log(0.0001),'priorSD',0.1,'NSamples',1,'type','real');

