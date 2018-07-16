lam1=[0.5 0.6 0.7 0.8 0.9];
acc1=[65.33 66.66 66.67 69.34 70];

lam2=[0.5 0.6 0.7 0.8 0.9 0.95];
acc2=[66 82.66 88 66 90.6 89.3];

lam3=[0.5 0.6 0.7 0.8 0.9 0.95];
acc3=[81.33 85.33 92.67 92.67 93.33 92.67];

lam4=[0.5 0.6 0.7 0.8 0.9 0.95];
acc4=[91.33 92 93.33 93.33 93.33 93.33];


clear y label;
box on
hold on 
graph1=plot(lam1,acc1,'green');
graph2=plot(lam2,acc2,'blue');
graph3=plot(lam3,acc3,'red');
graph4=plot(lam4,acc4,'black');
set(graph1,'LineWidth',1.5);
set(graph2,'LineWidth',1.5);
set(graph3,'LineWidth',1.5);
set(graph4,'LineWidth',1.5);
yticks([50 60 70 80 90 100])
scatter(lam1,acc1,'green');
scatter(lam2,acc2,'blue');
scatter(lam3,acc3,'red','filled');
scatter(lam4,acc4,'black','filled');
legend('Crisp','Type-1 AP','Interval Type-2 AP' , 'General Type-2 AP');
ylim([50 100]);
hold off