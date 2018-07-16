%x=[1 2 99 91 4000 4005]';
% % % delimiterIn =',';
 data=split(importdata('iris dataset.txt'),",");
 size(data)
% disp(data)
 X=data(:,1:4);
 out=data(:,5);
% disp(X)
% disp(Y)
 x=str2double(X);
%x= [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 52 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 79 80 81 82 83 84 85 86 88 89 90 91 92 93 94 95 96 97 98 99 100 102 107 114 120 122 127 128 139 143  51 53 77 78 87 101 103 104 105 106 108 109 110 111 112 113 115 116 117 118 119 121 123 124 125 126 129 130 131 132 133 134 135 136 137 138 140 141 142 144 145 146 147 148 149 150]';
%x=[3 4 3 2 1; 4 3 5 1 1;3 5 3 3 3; 2 1 3 3 2;1 1 3 2 3];
%x=[1 2 3 4 5 99 98 97 96]';
n=length(x);
p=size(x,2);
Z=x;
m=n;
mf=ones(n,p,n);
prevmf=ones(n,p,n);
dis=zeros(n,p,m);
disbar=zeros(n,m);
devbar=zeros(n,m);
dev=zeros(n,p,m);
beta=1;
prevdisbar=zeros(n,m);
comparison2=1;
s=zeros(n,m);

%initializing prevmf
for i=1:n
        prevmf(i,:,i)= 0;
        
end 
mf=prevmf;
%initializing dis and disbar
for i=1:n
        for k=1:m
            if(i~=k)
            sum2=0;
            sum3=0;
                for j=1:p
                    dis(i,j,k)=abs(x(i,j)-Z(k,j));
                    sum2 =sum2 + dis(i,j,k);
                end
          disbar(i,k)= sum2/p;
           end
        end
end

l=0;
%iterations
%while comparison2==1
for iter=1:50
    l=l+1;
    for i=1:n
        for k=1:m
            if(i~=k)
                sum2=0;
                sum3=0;
                    for j=1:p
                       %dis(i,j,k)=abs(x(i,j)-Z(k,j));
                        sum2 =sum2 + dis(i,j,k)*prevmf(i,j,k);
                        sum3= sum3+ prevmf(i,j,k);
                    end
                %if sum3~=0
                    disbar(i,k)= sum2/sum3;
                %else
                    %disbar(i,k)=0;
              
               % end
           end
        end
    end


%Dev and devbar
for i=1:n
  for k=1:m
      if(i~=k)
    sum4=0;
    for j=1:p
      dev(i,j,k)=abs(dis(i,j,k)-disbar(i,k));
      sum4=sum4+dev(i,j,k);
    end
    devbar(i,k)=sum4/p;
   end
   end
end


%membershipfunction
for i2=1:n
  for k=1:m
      if(i2~=k)
         for j=1:p
             if(devbar(i2,k)~=0)
             mf(i2,j,k)= exp(-dev(i2,j,k)^beta / devbar(i2,k)^beta);
             else
                 mf(12,j,k)=0;
             end
             
             %disp([dev(i2,j,k) devbar(i2,k) mf(i2,j,k)]);
         end
      end
   end
end


prevmf=mf;
comparison2=0;

for i=1:n
    for k=1:m
        comparison1=abs(prevdisbar(i,k)-disbar(i,k));
            if comparison1>0.01
                comparison2 =1;
            end 
    end
end
prevdisbar=disbar;


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
%
% I=(mean(mean(s))*eye(n))/4;
% s=s+I;
 for i=1:n
   s(i,i)= (mins - 2*(maxs-mins));
end
disp(s);
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
disp(idx)
% 
% for i=1:n
%     disp([i,idx(i)]);
% end

uidx=unique(idx);
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
disp(l);

disp(acc);

box on
color=(idx+44)/(sum(uidx));
scatter3((x(:,2)-min(x(:,2)))/(max(x(:,2))-min(x(:,2))),(x(:,3)-min(x(:,3)))/(max(x(:,3))-min(x(:,3))),(x(:,4)-min(x(:,4)))/(max(x(:,4))-min(x(:,4))),100,color,'filled');
