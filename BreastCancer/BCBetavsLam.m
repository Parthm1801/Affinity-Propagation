lam1=[0.9 0.91 0.92 0.93 0.94 0.95];
acc1=[65.24 72.1 85.69 82.69 85.26 65.52];

lam2=[0.9 0.91 0.92 0.93 0.94 0.95 ];
acc2=[86.83 85.56 85.69 86.83 86.83 86.83 ];

lam3=[0.86 0.87 0.88  0.89  0.9 0.91 0.92 0.93 0.94 0.95 ];
acc3=[94.84 91.41 88.13 89.13  90.4 86.83 86.56 85.98 86.26 85.54  ];

lam4=[0.86 0.87 0.88  0.89  0.9 0.91 0.92 0.93 0.94 0.95 ];
acc4=[91.56 91.56 91.56 89.27 89.27 89.12 89.12  89.41 89.12 89.12  ];

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
legend('Crisp','Type-1 AP','Interval Type-2 AP', 'General Type-2 AP');
ylim([50 100]);
hold off