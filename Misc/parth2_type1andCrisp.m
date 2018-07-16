%x=[1 2 99 91 4000 4005]';
% % % delimiterIn =',';
%  data=split(importdata('iris dataset.txt'),",");
%  size(data)
% % % disp(data)
%  X=data(:,1:4);
%  Y=data(:,5);
% % % disp(X)
% % % disp(Y)
%  x=str2double(X);
x=[3 4 3 2 1; 4 3 5 1 1;3 5 3 3 3; 2 1 3 3 2;1 1 3 2 3];
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
    for k=1:n
        if(i==k)
        prevmf(i,:,k)= 0;
        end
    end
end 

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


%iterations
while comparison2==1
    
    for i=1:n
        for k=1:m
            if(i~=k)
                sum2=0;
                sum3=0;
                    for j=1:p
                        dis(i,j,k)=abs(x(i,j)-Z(k,j));
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

disp(mf);

    

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

 for i=1:n
   s(i,i)= (mins - 0.1*(maxs-mins));
end
%disp(s);
S=s;

% Initialize messages
N=size(S,1); 
A=zeros(N,N); 
R=zeros(N,N); 
S=(S+(eps*randn(N,N))*(max(S(:))-min(S(:)))); % Remove degeneracies
lam=0.5; % Set damping factor
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

%disp([unique(idx),histc(idx,unique(idx))]);

