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
	disp(num2str(maxevals));
  for pind=1:numel(inproblems)
		inproblem = inproblems(pind);
		eval(['load ourresults/five-fixed-out-1-10-' inproblem.name '-' num2str(maxevals)]);
		disp(['problem: ' inproblem.name]);
		e = [];
		for i=1:numel(outproblems)
			e = [e; outproblems(i).allerrors([1:3 6])];
		end
		m = mean(e)
		s = std(e)
	end
end
