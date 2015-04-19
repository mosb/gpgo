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
%problem.name = 'griewank10';
%inproblems = [inproblems; problem];

%ackley 
problem.name = 'ackley2';
inproblems = [inproblems; problem];

%ackley 
problem.name = 'ackley5';
inproblems = [inproblems; problem];

%ackley 
%problem.name = 'ackley10';
%inproblems = [inproblems; problem];

%rast 2
problem.name = 'rast';
inproblems = [inproblems; problem];

inproblems = inproblems(probinds);

fid = fopen('egoetcresults3','w');  
for maxevals = [10]
  disp(['maxevals: ' num2str(maxevals)]);
  fprintf(fid,['maxevals: ' num2str(maxevals) '\n']);
  for pind=1:numel(inproblems)
    inproblem = inproblems(pind);
    eval(['load testfunctions/' inproblem.name]);
    disp(['problem: ' inproblem.name]);
    outproblems = [];
    for off = 1:10
      display(['offset number: ' num2str(off)]);
      problem = thisproblem(off);
      problem.maxevals = maxevals;
			
      Prob = glcAssign(problem.shortname,problem.lowrbnd,problem.upprbnd,problem.name,[],[],[],'',[],[],problem.x0);
      Prob.optParam.MaxFunc = maxevals * problem.n;
      Prob.optParam.IterPrint = 0;
      problem.Prob = Prob;
			
      results = ego(Prob);
      problem.egoresults = results;
      problem.egoerror = ...
          (problem.f(problem.x0) - results.f_k) / (problem.f(problem.x0) - problem.optimum);
      disp(['ego error: ' num2str(problem.egoerror)]);
			
      results = rbfSolve(Prob);
      problem.rbfresults = results;
      problem.rbferror = ...
          (problem.f(problem.x0) - results.f_k) / (problem.f(problem.x0) - problem.optimum);
      disp(['rbf error: ' num2str(problem.rbferror)]);
			
      results = glbSolve(Prob);
      problem.directresults = results;
      problem.directerror = ...
          (problem.f(problem.x0) - results.f_k) / (problem.f(problem.x0) - problem.optimum);
      disp(['direct error: ' num2str(problem.directerror)]);
			
      outproblems = [outproblems; problem];
      outstring = [problem.name ' ego: ' num2str(problem.egoerror) ...
                                ' rbf: ' num2str(problem.rbferror) ... 
                                ' direct: ' num2str(problem.directerror)];
      disp(outstring);
      fprintf(fid, [outstring '\n']);
    end
    cd tomlabresults
    eval(['save again-out-' problem.name '-' num2str(maxevals) ' outproblems']);
    cd '..'
  end
end
