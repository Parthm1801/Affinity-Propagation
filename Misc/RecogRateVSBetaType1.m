T1beta=[1 2 5 10 20 30 40 50 60 70 80 90 100];
T1recograte=[89.33 90 90 86 84 82.67 84 84 84 84.67 84.67 84 84.67];

%clear title, x label, y label;
hold on 
plot(T1beta,T1recograte,'red');
title('Type 1 : Recognition Rate vs Beta');
xlabel('Beta');
%ylabel('Recognition Rate(%)');
legend('Type-1 AP');
scatter(T1beta,T1recograte,'Red','filled');
ylim([0 100]);
hold off