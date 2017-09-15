% multivariate gaussian
function y = fitfun2(x, n)

%pause(1);

s= mvnpdf(x, [-5 -5]) + mvnpdf(x, [5 5]); 
y = log(s);

