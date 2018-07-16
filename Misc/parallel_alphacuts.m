data=split(importdata('iris dataset.txt'),",");
 size(data)
% disp(data)
 X=data(:,1:4);
 Y=data(:,5);
% disp(X)
% disp(Y)
 x=str2double(X);
%x=[1 2; 3 4;99 91; 98 99;4000 4005;4000 4002];
beta1=2;
beta2=5;
p=2:0.1:5;
c=(beta1+beta2)/2;
sig=(beta2-beta1)/7;
syms m;
gauss(m)=gaussmf(m,[sig c]);
boss=matlabFunction(gauss);
gaussinv=finverse(gauss);
F=matlabFunction(gaussinv);
L=zeros(5,2);
update=0.01;
for rownum=1:5
        L(rownum,1)=beta1+beta2-F(rownum/5);
        L(rownum,2)=F(rownum/5);
end

i=1;
j=1;
k=120;
mfup=1;
mfdown=0.99;

% for i=1:5
%     disp([s(i,1),s(i,2)]);
% end

% 
% hold on        
% plot(p,boss(p));
% for rownum=1:5
%     plot(L(rownum,1),boss(L(rownum,1)),'*');
%     plot(L(rownum,2),boss(L(rownum,2)),'*');
% end
% hold off

beta=L;
u=zeros(5,1);
for rownum=1:5
    u(rownum,1)=boss(F(rownum/5));
end

points=zeros(10,1);
alpha=zeros(10,1);
onemore=zeros(5,1);
ekaur=zeros(5,1);

parfor iter=1:5
%%%%%%%%%%%%%%%%%%%%%major loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=length(x);
p=size(x,2);
Z=x;
m=n;
    %Initalizations for membership function corresponding to beta1
mf1=ones(n,p,m);
prevmf1=ones(n,p,m);
dis1=zeros(n,p,m);
disbar1=zeros(n,m);
devbar1=zeros(n,m);
dev1=zeros(n,p,m);
prevdisbar1=zeros(n,m);
comparison2=1;


%Initializations for membership function corresponding to beta2
mf2=ones(n,p,m);
prevmf2=ones(n,p,m);
dis2=zeros(n,p,m);
disbar2=zeros(n,m);
devbar2=zeros(n,m);
dev2=zeros(n,p,m);
prevdisbar2=zeros(n,m);
comparison4=1;


% %Initializations for overall membership functions
% dis=zeros(n,p,m);
% disbar=zeros(n,m);
% s=zeros(n,m);


%Computations for Membership function corresponding to beta1

%initializing prevmf1
for rownum=1:n
    for featurenum=1:j
        
        prevmf1(rownum,featurenum,rownum)= 0;
        
    end
end 

%initializing dis1 and disbar1
for rownum=1:n
        for clusternum=1:m
            if(rownum~=clusternum)
            sum2=0;
            sum3=0;
                for featurenum=1:p
                    dis1(rownum,featurenum,clusternum)=abs(x(rownum,featurenum)-Z(clusternum,featurenum));
                    sum2 =sum2 + dis1(rownum,featurenum,clusternum);
                end
          disbar1(rownum,clusternum)= sum2/p;
           end
        end
end

dis=dis1;
%iterations
for l=1:20
    
    for rownum=1:n
        for clusternum=1:m
            if(rownum~=clusternum)
                sum2=0;
                sum3=0;
                    for featurenum=1:p
                        dis1(rownum,featurenum,clusternum)=abs(x(rownum,featurenum)-Z(clusternum,featurenum));
                        sum2 =sum2 + dis1(rownum,featurenum,clusternum)*prevmf1(rownum,featurenum,clusternum);
                        sum3= sum3+ prevmf1(rownum,featurenum,clusternum);
                    end
                if sum3~=0
                    disbar1(rownum,clusternum)= sum2/sum3;
                else
                    disbar1(rownum,clusternum)=0;
              
                end
           end
        end
    end


%Dev1 and devbar1
for rownum=1:n
  for clusternum=1:m
      if(rownum~=clusternum)
    sum4=0;
    for featurenum=1:p
      dev1(rownum,featurenum,clusternum)=abs(dis1(rownum,featurenum,clusternum)-disbar1(rownum,clusternum));
      sum4=sum4+dev1(rownum,featurenum,clusternum);
    end
    devbar1(rownum,clusternum)=sum4/p;
   end
   end
end


%membershipfunction1
for i2=1:n
  for clusternum=1:m
      if(i2~=clusternum)
         for featurenum=1:p
             if(devbar1(i2,clusternum)~=0)
             mf1(i2,featurenum,clusternum)= exp(-dev1(i2,featurenum,clusternum)^beta(iter,1) / devbar1(i2,clusternum)^beta(iter,1));
             else
                 mf1(i2,featurenum,clusternum)=0;
             end
         end
      end
   end
end


prevmf1=mf1;
comparison2=0;

for rownum=1:n
    for clusternum=1:m
        comparison1=abs(prevdisbar1(rownum,clusternum)-disbar1(rownum,clusternum));
            if comparison1>0.01
                comparison2 =1;
            end 
    end
end
prevdisbar1=disbar1;

end





%Computations for Membership function corresponding to beta2

%initializing prevmf2
for rownum=1:n
    for clusternum=1:n
        if(rownum==clusternum)
        prevmf2(rownum,:,clusternum)= 0;
        end
    end
end 

%initializing dis2 and disbar2
for rownum=1:n
        for clusternum=1:m
            if(rownum~=clusternum)
            sum2=0;
            sum3=0;
                for featurenum=1:p
                    dis2(rownum,featurenum,clusternum)=abs(x(rownum,featurenum)-Z(clusternum,featurenum));
                    sum2 =sum2 + dis2(rownum,featurenum,clusternum);
                end
          disbar2(rownum,clusternum)= sum2/p;
           end
        end
end


%iterations
%while comparison2==1
 for l=1:20   
    for rownum=1:n
        for clusternum=1:m
            if(rownum~=clusternum)
                sum2=0;
                sum3=0;
                    for featurenum=1:p
                        dis2(rownum,featurenum,clusternum)=abs(x(rownum,featurenum)-Z(clusternum,featurenum));
                        sum2 =sum2 + dis2(rownum,featurenum,clusternum)*prevmf2(rownum,featurenum,clusternum);
                        sum3= sum3+ prevmf2(rownum,featurenum,clusternum);
                    end
                if sum3~=0
                    disbar2(rownum,clusternum)= sum2/sum3;
                else
                    disbar2(rownum,clusternum)=0;
              
                end
           end
        end
    end


%Dev2 and devbar2
for rownum=1:n
  for clusternum=1:m
      if(rownum~=clusternum)
    sum4=0;
    for featurenum=1:p
      dev2(rownum,featurenum,clusternum)=abs(dis2(rownum,featurenum,clusternum)-disbar2(rownum,clusternum));
      sum4=sum4+dev2(rownum,featurenum,clusternum);
    end
    devbar2(rownum,clusternum)=sum4/p;
   end
   end
end


%membershipfunction2
for i2=1:n
  for clusternum=1:m
      if(i2~=clusternum)
         for featurenum=1:p
             if(devbar2(i2,clusternum)~=0)
             mf2(i2,featurenum,clusternum)= exp(-dev2(i2,featurenum,clusternum)^beta(iter,2) / devbar2(i2,clusternum)^beta(iter,2));
             else
                 mf2(i2,featurenum,clusternum)=0;
             end
         end
      end
   end
end


prevmf2=mf2;
comparison4=0;

for rownum=1:n
    for clusternum=1:m
        comparison3=abs(prevdisbar2(rownum,clusternum)-disbar2(rownum,clusternum));
            if comparison3>0.01
                comparison4 =1;
            end 
    end
end
prevdisbar2=disbar2;

 end
 

     points(iter,1)= mf1(i,j,k);
     onemore(iter,1)= mf2(i,j,k);
     alpha(iter,1)=u(iter);
     ekaur(iter,1)=u(iter);
 
 %%%%%%%%%%%%%%%%%%%major loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
end

for iter=1:5
    points(11-iter,1)=onemore(iter,1);
    alpha(11-iter,1)=ekaur(iter,1);
end
  points(5)=[];
  alpha(5)=[];

count=0;
var=points(1);
for i=1:9
    if points(i)==var
        count=count+1;
    end
end

if count==9
    secgrade=0;  %%%%%%%%% but actually it should be 1 at x axis =points(i)
else
secgrade = fit(points,alpha,'gauss3');
end

mftemp=mfdown;

if count~=9
     while mfdown<=mftemp && mftemp<=mfup
          disp("hey");
        sum5= sum5+ mftemp*secgrade(mftemp);
        sum6= sum6 + secgrade(mftemp);
        mftemp=mftemp+update;
      end
else
    sum5=1;
    sum6=1;
end
disp("yo")



