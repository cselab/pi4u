% multivariate gaussian
function y = fitfun1(x, n)

%pause(1);

s= mvnpdf(x); 
y = log(s);

