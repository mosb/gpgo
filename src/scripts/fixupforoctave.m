format compact
warning off
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

%inproblems = inproblems(probinds);
maxevals = 10;

for pind=1:numel(inproblems)
	inproblem = inproblems(pind);
	disp(inproblem.name)
	eval(['load out-' inproblem.name '-' num2str(maxevals)]);
	eval(['save out-' inproblem.name '-' num2str(maxevals) ' -v6']);
end
