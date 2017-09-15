function display_gen(datafile,Nth,i,j)

c=load(datafile); 
scatter(c(:,i),c(:,j),10,c(:,Nth+1))
xlabel('var1')
ylabel('var2')
colorbar

end

	


