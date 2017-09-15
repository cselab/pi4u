% data-driven inference
function y = fitfun3(x, n)

a = x(1);
sigma = x(2);

data = load('data3.txt');
N=size(data,1);

fx = data(:,1);
fy = data(:,2);
	
y = a*fx;

SSE = sumsqr(y-fy);
sy=sigma;
logn = -0.5*N*log(2*pi)-0.5*N*log(sy*sy)-0.5*SSE/(sy*sy);

y=logn;
