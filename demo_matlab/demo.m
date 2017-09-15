function demo()

path(path,'../lib_matlab')

% input arguments
Nth = 2
MaxStages = 20
PopSize = 4096
id = 2		% unique identifier, appended to curgen_db files
lb = [-10 -10]
%lb = [0 0]
ub = [10 10]

logEv = tmcmc('fitfun3',id,Nth,MaxStages,PopSize,lb,ub)

% addpath('../lib_matlab');
% display_gen('curgen_db_002_008.txt',2,1,2)

