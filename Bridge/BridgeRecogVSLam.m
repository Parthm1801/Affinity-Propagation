lam1=[ 0.6 0.7 0.8 0.9  ];
acc1=[99.14 99.14 93.96 78.44];

lam2=[0.6 0.7 0.8 0.9];
acc2=[99.14 99.14 99.14 94.82];

clear y label;
box on
hold on 
graph1=plot(lam1,acc1,'blue');
graph2=plot(lam2,acc2,'red');
set(graph1,'LineWidth',1.5);
set(graph2,'LineWidth',1.5);
yticks([50 60 70 80 90 100])
scatter(lam1,acc1,'blue');
scatter(lam2,acc2,'red', 'filled');
legend('Type-1 AP','Interval Type-2 AP');
ylim([50 100]);
hold off