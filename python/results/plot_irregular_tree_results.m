load('regular_tree_results_1000trials.mat')

[num_d,max_T] = size(pd_ml);

figure(1)
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


figure(3)
[m,I] = min(pd_ml);
plot(1:max_T,d_values(I)+1)
xlabel('Timestep (T)')
ylabel('Optimal d_o')


