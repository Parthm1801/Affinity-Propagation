betaIT2=[1 10 20 30 40 50 60 70 80 90 100];
acc2=[93.33 88.33 87.33 85.33 86 87.33 90 86 86.67 86.67 86.67];

betaT1=[1 10 20 30 40 50 60 70 80 90 100];
acc1=[89.33  86 84 82.67 84 84 84 84.67 84.67 84 84.67];

betaGT2=[1 10 20 30 40 50 60 70 80 90 100];
acc3=[93.33 89.33 87.33 89.33 88.67 88.67  89.33  88.67 88.67 88.67 88.67 ];


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
legend('Type-1 AP','Interval Type-2 AP' , 'General Type-2 AP');
ylim([50 100]);

hold off