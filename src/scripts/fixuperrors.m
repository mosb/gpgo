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
problem.periodic = true;
inproblems = [inproblems; problem];
%griewank
problem.name = 'griewank2';
problem.periodic = true;
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
problem.periodic = true;
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
problem.periodic = true;
inproblems = [inproblems; problem];

inproblems = inproblems(probinds);

for maxevals = [10]
    disp(['maxevals: ' num2str(maxevals)]);
    for pind=1:numel(inproblems)
        inproblem = inproblems(pind);
        eval(['load ourresults/five-out-1-10-' inproblem.name '-' num2str(maxevals)]);
        oldproblems = outproblems([1:min(numel(outproblems),10)]);
        disp(num2str(numel(oldproblems)));
        clear outproblems 
        
        for i=1:numel(oldproblems)
            disp(['problem: ' inproblem.name]);
            problem = oldproblems(i);
            problem.egoerror = ...
                (problem.forig(problem.x0) - problem.egoresults.f_k) / ...
                (problem.forig(problem.x0) - problem.optimum);
            problem.rbferror = ...
                (problem.forig(problem.x0) - problem.rbfresults.f_k) / ...
                (problem.forig(problem.x0) - problem.optimum);
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
            problem.allerrors = [problem.egoerror problem.rbferror ...
                                problem.directerror problem.ourerrors'];
            outproblems(i) = problem;
            
            cd ourresults
            eval(['save five-fixed-out-1-10-' problem.name '-' num2str(maxevals) ' outproblems off']);
            cd '..'
        end
    end
end
