betaT1=[10 20 30 40 50 60 70 80 90 100];
acc1=[86.12 86.84 85.85 85.84 85.56 85.98 85.84 85.98 85.98 85.98];

betaIT2=[10 20 30 40 50 60 70 80 90 100];
acc2=[90.41 91.41 89.99 90.13 90.7 89.98 90.13 90.56 90.56 90.13];

betaGT2=[10 20 30 40 50 60 70 80 90 100];
acc3=[91.41 91.41 91.41 91.41 91.41 91.41 91.41 91.41 91.41 91.56];

clear y label;
box  on
hold on 

graph1=plot(betaT1,acc1,'blue');
graph2=plot(betaIT2,acc2,'red');
graph3=plot(betaGT2,acc3,'black');
set(graph1,'LineWidth',1.5);
set(graph2,'LineWidth',1.5);
set(graph3,'LineWidth',1.5);
yticks([50 60 70 80 90 100])
scatter(betaT1,acc1,30,'blue');
scatter(betaIT2,acc2,30,'red','filled');
scatter(betaGT2,acc3,30,'black','filled');
legend('Type-1 AP','Interval Type-2 AP', 'General Type-2 AP');
%legend('Type-2 AP','Type-2 Data','Type-1 AP','Type-1 Data');
ylim([50 100]);

hold off
