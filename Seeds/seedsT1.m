data=importdata('seeds.txt');
 x=data(:,1:7);
 out=data(:,8);
 
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
beta=0.5;
prevdisbar=zeros(n,m);
comparison2=1;
s=zeros(n,m);

%initializing prevmf
for i=1:n
        prevmf(i,:,i)= 0;
  
end

mf=prevmf;
%initializing dis and disbar
for J=1:n
    dis(:,:,J)=abs(x-Z(J,:));
end
        disbar=squeeze(sum(dis,2)/p);


%iterations
for it=1:50
    
sum2=sum((dis.*mf),2);
sum3=sum(mf,2);
disbar=sum2./sum3;


%Dev and devbar
B=repmat(disbar,1,p,1);
dev=abs(dis-B);
devbar=(sum(dev,2))/p;

P=repmat(devbar,1,p,1);

%membershipfunction1
prevmf=exp(1).^(-(dev.^beta)./(P.^beta));
prevmf(isnan(prevmf))=0;
             


mf=prevmf;
% comparison2=0;
% 
% G=abs(prevdisbar-disbar);

% if sum(G(:)>0.01)~=0
%     comparison2 =1;
% end

prevdisbar=disbar;


end

    

%Similarity matrix
s=-squeeze(disbar);

for i=1:n
    s(i,i)=0;
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
   s(i,i)= (mins - 3.7*(maxs-mins));
end

S=s;

% Initialize messages
N=size(S,1); 
A=zeros(N,N); 
R=zeros(N,N); 
S=(S+(eps*randn(N,N))*(max(S(:))-min(S(:))))*40; % Remove degeneracies
lam=0.95; % Set damping factor
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
% disp(idx)

for i=1:n
    disp([i,idx(i)]);
end

disp([unique(idx),histc(idx,unique(idx))]);
acc=zeros(3);
y=zeros(1,n);
ybar=zeros(1,n);

for i=1:n
    if (idx(i)>=1 && idx(i)<=70)
        ybar(i)=1;
    elseif (idx(i)>=71 && idx(i)<=140)
        ybar(i)=2;
    elseif (idx(i)>=141 && idx(i)<=210)
        ybar(i)=3;
    end
end
for i=1:n
    acc(out(i),ybar(i))= acc(out(i),ybar(i))+1;
end

disp(acc);

% color=(idx+20)/(sum(uidx)+200);
% %color=[color color color]+[rand() rand() rand()]
% scatter(x(:,1),x(:,2),100,color,'filled');
% set(gca,'Color','k');