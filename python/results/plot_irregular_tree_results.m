load('irregular_tree_results.mat')

[num_d,max_T] = size(pd_ml);

plot(num_infected', pd_ml')
xlabel('Number of nodes')
ylabel('Pd over 1000 trials')

figure(2)
integrated = []
for i = 1:num_d
    integrated = [integrated; cumtrapz(num_infected(i,:),pd_ml(i,:))];
end

plot(d_values+1,integrated(:,end))
xlabel('Assumed Degree d_o')
ylabel('Approximate Integral of Pd over N')