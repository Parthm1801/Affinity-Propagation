data=split(importdata('iris dataset.txt'),",");
 size(data)
% disp(data)
 X=data(:,1:4);
 out=data(:,5);
% disp(X)
% disp(Y)
 x=str2double(X);
n=length(x);
p=size(x,2);
Z=x;
m=n;


%Initalizations for membership function corresponding to beta1
mf1=ones(n,p,n);
prevmf1=ones(n,p,n);
dis1=zeros(n,p,m);
disbar1=zeros(n,m);
devbar1=zeros(n,m);
dev1=zeros(n,p,m);
beta1=1;
beta2=2.5;
prevdisbar1=zeros(n,m);
comparison2=1;


%Initializations for membership function corresponding to beta2
mf2=ones(n,p,n);
prevmf2=ones(n,p,n);
dis2=zeros(n,p,m);
disbar2=zeros(n,m);
devbar2=zeros(n,m);
dev2=zeros(n,p,m);

prevdisbar2=zeros(n,m);
comparison4=1;


%Initializations for overall membership functions
dis=zeros(n,p,m);
disbar=zeros(n,m);
s=zeros(n,m);
mfdown=zeros(n,p,m);
mfup=zeros(n,p,m);
mftemp=zeros(n,p,m);
update=0.01;


%Computations for Membership function corresponding to beta1

%initializing prevmf1
for i=1:n
       prevmf1(i,:,i)= 0;
end 

%initializing dis1 and disbar1
for J=1:n
    dis1(:,:,J)=abs(x-Z(J,:));
end
        disbar1=squeeze(sum(dis1,2)/p);
dis=dis1;
%iterations
 %while comparison2==1
for l=1:29
   sum21=sum((dis1.*prevmf1),2);
sum31=sum(prevmf1,2);
sum31=(sum31==0)+sum32;
disbar1=sum21./sum31;
disbar1(isnan(disbar1))=0;

sum22=sum((dis1.*prevmf2),2);
sum32=sum(prevmf2,2);
sum32=(sum32==0)+sum32;
disbar2=sum22./sum32;
disbar2(isnan(disbar2))=0;

B1=repmat(disbar1,1,p,1);
B2=repmat(disbar2,1,p,1);

dev1=abs(dis1-B1);
dev2=abs(dis1-B2);

devbar1=(sum(dev1,2))/p;
devbar2=(sum(dev2,2))/p;

P1=repmat(devbar1,1,p,1);
P2=repmat(devbar2,1,p,1);

    mf1=exp(1).^(-(dev1.^beta1)./P1.^beta1);
    mf2=exp(1).^(-(dev2.^beta2)./P2.^beta2);
    mf1(isnan(mf1))=0;
    mf2(isnan(mf2))=0;
    
prevmf1=mf1;
prevmf2=mf2;
comparison2=0;

G1=abs(prevdisbar1-disbar1);
G2=abs(prevdisbar2-disbar2);

if sum(sum(sum(G1>0.01)))~=0 || sum(sum(sum(G1>0.01)))~=0
    comparison2=1;
end

prevdisbar1=disbar1;
prevdisbar2=disbar2;

end


%Computing overall membership function

mfdown=(mf1<=mf2).*mf1 + (mf1>mf2).*mf2;
mfup=(mf1<=mf2).*mf2 + (mf1>mf2).*mf1;

mf=SG(x,mfup,mfdown);



%Computing dis and disbar
disbar= sum((dis.*mf),2)./sum(mf,2);
               


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



uidx=unique(idx);
for i=1:length(uidx)
    idx(uidx(i))=uidx(i);
end
for i=1:n
    disp([i,idx(i)]);
end
uidx=unique(idx);
for i=1:length(uidx)
    idx(uidx(i))=uidx(i);
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

disp(l);
color=(idx+44)/(sum(uidx));

scatter3((x(:,2)-min(x(:,2)))/(max(x(:,2))-min(x(:,2))),(x(:,3)-min(x(:,3)))/(max(x(:,3))-min(x(:,3))),(x(:,4)-min(x(:,4)))/(max(x(:,4))-min(x(:,4))),100,color,'filled');

%%%%%%%%%%% secondary membership fn %%%%%%%%%%%%%%%%%%


function mf=SG(x,mfup,mfdown)
beta1=1;
beta2=2.5;
p=beta1:0.1:beta2;
c=(beta1+beta2)/2;
sig=(beta2-beta1)/7;
syms m;
gauss(m)=gaussmf(m,[sig c]);
boss=matlabFunction(gauss);
gaussinv=finverse(gauss);
F=matlabFunction(gaussinv);
L=zeros(3,2);
update=0.01;
for I=1:3
        L(I,1)=beta1+beta2-F(I/4);
        L(I,2)=F(I/4);
end
n=length(x);
p=size(x,2);
Z=x;
m=n;


beta=L;
u=zeros(3,1);
for I=1:3
    u(I,1)=boss(F(I/4));
end

points=zeros(6,n,p,n);
alpha=zeros(6,n,p,n);
for iter=1:3
%%%%%%%%%%%%   major      %%%%%%%%%%%%%%%%%%%%%%%

    %Initalizations for membership function corresponding to beta1
mf1=ones(n,p,m);
prevmf1=ones(n,p,m);
dis1=zeros(n,p,m);
disbar1=zeros(n,p,m);
devbar1=zeros(n,m);
dev1=zeros(n,p,m);
prevdisbar1=zeros(n,m);
comparison2=1;

%Initializations for membership function corresponding to beta2
mf2=ones(n,p,m);
prevmf2=ones(n,p,m);

disbar2=zeros(n,p,m);
devbar2=zeros(n,m);
dev2=zeros(n,p,m);

comparison4=1;

%initializing dis and disbar

for I=1:n
    prevmf1(I,:,I)=0;
    prevmf2(I,:,I)=0;
end
for J=1:n
    dis1(:,:,J)=abs(x-Z(J,:));
end
        prevdisbar1=squeeze(sum(dis1,2)/p);
        prevdisbar2=prevdisbar1;

%iterations
for l=1:3
    
sum21=sum((dis1.*prevmf1),2);
sum31=sum(prevmf1,2);
sum31=(sum31==0)+sum31;
disbar1=sum21./sum31;
disbar1(isnan(disbar1))=0;

sum22=sum((dis1.*prevmf2),2);
sum32=sum(prevmf2,2);
sum32=(sum32==0)+sum32;
disbar2=sum22./sum32;
disbar2(isnan(disbar2))=0;

B1=repmat(disbar1,1,p,1);
B2=repmat(disbar2,1,p,1);

dev1=abs(dis1-B1);
dev2=abs(dis1-B2);

devbar1=(sum(dev1,2))/p;
devbar2=(sum(dev2,2))/p;

P1=repmat(devbar1,1,p,1);
P2=repmat(devbar2,1,p,1);

    mf1=exp(1).^(-(dev1.^beta(iter,1))./P1.^beta(iter,1));
    mf2=exp(1).^(-(dev2.^beta(iter,2))./P2.^beta(iter,2));
    mf1(isnan(mf1))=0;
    mf2(isnan(mf2))=0;


prevmf1=mf1;
prevmf2=mf2;
comparison2=0;

% G1=abs(prevdisbar1-disbar1);
% G2=abs(prevdisbar2-disbar2);
% 
% if sum(sum(sum(G1>0.01)))~=0 || sum(sum(sum(G2>0.01)))~=0
%     comparison2=1;
% end
% 
prevdisbar1=disbar1;
prevdisbar2=disbar2;

end
     points(l,:,:,:)= mf1;
     points(7-l,:,:,:)= mf2;
     
    alpha(l,:,:,:)=repmat(u(l),n,p,n);
    alpha(7-l,:,:,:)=repmat(u(l),n,p,n);
end

slope1=zeros(n,p,n);
slope2=zeros(n,p,n);
slope3=zeros(n,p,n);
slope4=zeros(n,p,n);

slope1=(alpha(2,:,:,:)-alpha(1,:,:,:))./(points(2,:,:,:)-points(1,:,:,:));
slope2=(alpha(3,:,:,:)-alpha(2,:,:,:))./(points(3,:,:,:)-points(2,:,:,:));
slope3=(alpha(4,:,:,:)-alpha(5,:,:,:))./(points(4,:,:,:)-points(5,:,:,:));
slope4=(alpha(5,:,:,:)-alpha(6,:,:,:))./(points(5,:,:,:)-points(6,:,:,:));


mftemp=mfdown;
sum5=zeros(n,p,n);
sum6=sum5;

for i=1:n
    tic
    for k=1:n
        if i~=k
            for j=1:p
                while mfdown(i,j,k)<=mftemp(i,j,k) && mftemp(i,j,k)<=mfup(i,j,k)
                 line1=slope1(1,i,j,k)*(mftemp(i,j,k)-points(1,i,j,k))+alpha(1,i,j,k);
                 line2=slope2(1,i,j,k)*(mftemp(i,j,k)-points(3,i,j,k))+alpha(3,i,j,k);
                 line3=slope3(1,i,j,k)*(mftemp(i,j,k)-points(4,i,j,k))+alpha(4,i,j,k);
                 line4=slope4(1,i,j,k)*(mftemp(i,j,k)-points(6,i,j,k))+alpha(6,i,j,k);
                 mint=min(min(min(line1,line2),min(line3,line4)),1);
                 
                 sum5(i,j,k)= sum5(i,j,k)+ mftemp(i,j,k)*max(0,mint);
                 sum6(i,j,k)= sum6(i,j,k) + max(0,mint);
                 mftemp(i,j,k)=mftemp(i,j,k)+update;
                
                end
            end
        end
    end
    toc
end

sum6=(sum6<0.0001)+sum6;
mf=sum5./(sum6.*(mfup-mfdown));
mf(isnan(mf))=0;
end