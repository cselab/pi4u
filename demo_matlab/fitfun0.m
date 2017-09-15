function y = fitfun0(x, n)

%pause(1);
sum = 0;
for j = 1:n-1;
    sum = sum+100*(x(j)^2-x(j+1))^2+(x(j)-1)^2;
end
y = -sum;

