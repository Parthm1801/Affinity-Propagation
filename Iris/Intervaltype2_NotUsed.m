%x=[1 2 3 4 99 97 98 ]';
 data=split(importdata('iris dataset.txt'),",");
 size(data)
% disp(data)
 X=data(:,1:4);
 Y=data(:,5);
% disp(X)
% disp(Y)
 x=str2double(X);
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
beta1=1;
prevdisbar1=zeros(n,m);
comparison2=1;


%Initializations for membership function corresponding to beta2
mf2=ones(n,p,m);
prevmf2=ones(n,p,m);
dis2=zeros(n,p,m);
disbar2=zeros(n,m);
devbar2=zeros(n,m);
dev2=zeros(n,p,m);
beta2=2.5;
prevdisbar2=zeros(n,m);
comparison4=1;


%Initializations for overall membership functions
dis=zeros(n,p,m);
disbar=zeros(n,m);
s=zeros(n,m);


%Computations for Membership function corresponding to beta1

%initializing prevmf1
for i=1:n
    for k=1:n
        if(i==k)
        prevmf1(i,:,k)= 0;
        end
    end
end 

%initializing dis1 and disbar1
for i=1:n
        for k=1:m
            if(i~=k)
            sum2=0;
            sum3=0;
                for j=1:p
                    dis1(i,j,k)=abs(x(i,j)-Z(k,j));
                    sum2 =sum2 + dis1(i,j,k);
                end
          disbar1(i,k)= sum2/p;
           end
        end
end

dis=dis1;
%iterations
for l=1:50
    
    for i=1:n
        for k=1:m
            if(i~=k)
                sum2=0;
                sum3=0;
                    for j=1:p
                        dis1(i,j,k)=abs(x(i,j)-Z(k,j));
                        sum2 =sum2 + dis1(i,j,k)*prevmf1(i,j,k);
                        sum3= sum3+ prevmf1(i,j,k);
                    end
                if sum3~=0
                    disbar1(i,k)= sum2/sum3;
                else
                    disbar1(i,k)=0;
              
                end
           end
        end
    end


%Dev1 and devbar1
for i=1:n
  for k=1:m
      if(i~=k)
    sum4=0;
    for j=1:p
      dev1(i,j,k)=abs(dis1(i,j,k)-disbar1(i,k));
      sum4=sum4+dev1(i,j,k);
    end
    devbar1(i,k)=sum4/p;
   end
   end
end


%membershipfunction1
for i2=1:n
  for k=1:m
      if(i2~=k)
         for j=1:p
             if(devbar1(i2,k)~=0)
             mf1(i2,j,k)= exp(-dev1(i2,j,k)^beta1 / devbar1(i2,k)^beta1);
             else
                 mf1(i2,j,k)=0;
             end
         end
      end
   end
end


prevmf1=mf1;
comparison2=0;

for i=1:n
    for k=1:m
        comparison1=abs(prevdisbar1(i,k)-disbar1(i,k));
            if comparison1>0.01
                comparison2 =1;
            end 
    end
end
prevdisbar1=disbar1;

end





%Computations for Membership function corresponding to beta2

%initializing prevmf2
for i=1:n
    for k=1:n
        if(i==k)
        prevmf2(i,:,k)= 0;
        end
    end
end 

%initializing dis2 and disbar2
for i=1:n
        for k=1:m
            if(i~=k)
            sum2=0;
            sum3=0;
                for j=1:p
                    dis2(i,j,k)=abs(x(i,j)-Z(k,j));
                    sum2 =sum2 + dis2(i,j,k);
                end
          disbar2(i,k)= sum2/p;
           end
        end
end


%iterations
%while comparison2==1
 for l=1:50   
    for i=1:n
        for k=1:m
            if(i~=k)
                sum2=0;
                sum3=0;
                    for j=1:p
                        dis2(i,j,k)=abs(x(i,j)-Z(k,j));
                        sum2 =sum2 + dis2(i,j,k)*prevmf2(i,j,k);
                        sum3= sum3+ prevmf2(i,j,k);
                    end
                if sum3~=0
                    disbar2(i,k)= sum2/sum3;
                else
                    disbar2(i,k)=0;
              
                end
           end
        end
    end


%Dev2 and devbar2
for i=1:n
  for k=1:m
      if(i~=k)
    sum4=0;
    for j=1:p
      dev2(i,j,k)=abs(dis2(i,j,k)-disbar2(i,k));
      sum4=sum4+dev2(i,j,k);
    end
    devbar2(i,k)=sum4/p;
   end
   end
end


%membershipfunction2
for i2=1:n
  for k=1:m
      if(i2~=k)
         for j=1:p
             if(devbar2(i2,k)~=0)
             mf2(i2,j,k)= exp(-dev2(i2,j,k)^beta2 / devbar2(i2,k)^beta2);
             else
                 mf2(i2,j,k)=0;
             end
         end
      end
   end
end


prevmf2=mf2;
comparison4=0;

for i=1:n
    for k=1:m
        comparison3=abs(prevdisbar2(i,k)-disbar2(i,k));
            if comparison3>0.01
                comparison4 =1;
            end 
    end
end
prevdisbar2=disbar2;

end

disp(mf1);
disp(mf2);

%Computing overall membership function
mf=(mf1+mf2)*0.5;

%Computing dis and disbar
for i=1:n
   for k=1:m
       if(i~=k)
            sum2=0;
            sum3=0;
                for j=1:p
                    sum2 =sum2 + dis(i,j,k)*mf(i,j,k);
                    sum3= sum3+ mf(i,j,k);
                end
                if sum3~=0
                    disbar(i,k)= sum2/sum3;
                else
                    disbar(i,k)=0;
              
                end
        end
    end
end


%Similarity matrix
for i=1:n
    for k=1:m
        if(i~=k)
            s(i,k)=-disbar(i,k);
        end
    end
end


%calculating max and min of similarity matrix

mins=min(min(s));
maxs=s(1,2);
for i=1:n
    for k=1:m
        if(i~=k && s(i,k)>=maxs)
            maxs=s(i,k);
        end
    end
end

            
%preference values
 for i=1:n
   s(i,i)= (mins - 2*(maxs-mins));
end

S=s;


% Initialize messages
N=size(S,1); 
A=zeros(N,N); 
R=zeros(N,N); 
S=(S+(eps*randn(N,N))*(max(S(:))-min(S(:)))); % Remove degeneracies
lam=0.9; % Set damping factor
for i=1:100
    % Compute responsibilities
    Rold=R;
    AS=A+S; 
    [Y,I]=max(AS,[],2);
    for i=1:N 
      AS(i,I(i))=-realmax; 
    end
    [Y2,I2]=max(AS,[],2);
    R=S-repmat(Y,[1,N]);
    for i=1:N 
      R(i,I(i))=S(i,I(i))-Y2(i);
    end
    R=(1-lam)*R+lam*Rold; % Dampen responsibilities

    % Compute availabilities
    Aold=A;
    Rp=max(R,0); 
    for k=1:N 
      Rp(k,k)=R(k,k); 
    end
    A=repmat(sum(Rp,1),[N,1])-Rp;
    dA=diag(A); 
    A=min(A,0); 
    for k=1:N 
      A(k,k)=dA(k);
    end
    A=(1-lam)*A+lam*Aold; % Dampen availabilities
end
E=R+A; % Pseudomarginals
I=find(diag(E)>0); 
K=length(I); % Indices of exemplars
[tmp,c]=max(S(:,I),[],2); 
idx=I(c); % Assignments
%disp(idx)

for i=1:n
    disp([i,idx(i)]);
end

disp([unique(idx),histc(idx,unique(idx))]);
acc=zeros(3,3);
y=zeros(1,n);
ybar=zeros(1,n);
for i=1:n
    if(out(i)== "Iris-setosa")
        y(i)=1;
    elseif (out(i)=="Iris-versicolor")
        y(i)=2;
    else
        y(i)=3;
    end
end
for i=1:n
    if (idx(i)>=1 && idx(i)<=50)
        ybar(i)=1;
    elseif (idx(i)>=51 && idx(i)<=100)
        ybar(i)=2;
    elseif (idx(i)>=101 && idx(i)<=150)
        ybar(i)=3;
    end
end
for i=1:n
    acc(y(i),ybar(i))= acc(y(i),ybar(i))+1;
end

disp(acc);

% clust=unique(idx);
% % colors=rand(clust(size,1),3);
% size(x)
% % plot(x(1),idx(1),'*','color',[idx(1) id)
% colors=[]
% for i=1:size(clust,1)
%     colors=[colors;repmat([rand() rand() rand()],sum(idx==clust(i)),1)]
% end
% scatter(x,idx,10,colors,'filled')