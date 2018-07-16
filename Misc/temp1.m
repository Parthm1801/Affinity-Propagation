x=[1 2 3; 56 78 65; 34 5 4; 69 31 40; 3 6 8; 7 19 23; 29 18 24; 20 100 4; 54 72 78 ; 60 98 30];
%x=randi(100,10,2);
%x=gpuArray(x);
n=length(x);
p=size(x,2);
Z=x;
m=n;
mf=ones(n,p,n);
%mf=gpuArray(mf);
prevmf=ones(n,p,n);
%prevmf=gpuArray(prevmf);
dis=zeros(n,p,m);
%dis=gpuArray(dis);
disbar=zeros(n,m);
%disbar=gpuArray(disbar);
devbar=zeros(n,m);
%devbar=gpuArray(devbar);
dev=zeros(n,p,m);
%dev=gpuArray(dev);
beta=1;
prevdisbar=zeros(n,m);
%prevdisbar=gpuArray(prevdisbar);
comparison2=1;
s=zeros(n,m);
%s=gpuArray(s);
% max=0;

disp("Hey");
while comparison2==1

%iterations
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
          disbar(i,k)= sum2/sum3;
           end
        end
    end
%  disp(dev);
%  disp(disbar);   
%     
%Dev and devbar
for i1=1:n
  for k=1:m
      if(i1~=k)
    sum4=0;
    for j=1:p
      dev(i1,j,k)=abs(dis(i1,j,k)-disbar(i1,k));
      sum4=sum4+dev(i1,j,k);
    end
    devbar(i1,k)=sum4/p;
   end
   end
end


%disp(devbar)

%membershipfunction
for i2=1:n
  for k=1:m
      if(i2~=k)
         for j=1:p
             mf(i2,j,k)= exp(-dev(i2,j,k)^beta / devbar(i2,k)^beta);
             %disp([dev(i2,j,k) devbar(i2,k) mf(i2,j,k)]);
         end
      end
   end
end

disp("**************************************************************************");
 % disp(prevmf)
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

for i=1:n
    s(i,i)= mins - m*(maxs-mins);
end

disp(s);

preva=zeros(n,n);
prevr=zeros(n,n);
a=zeros(n,n);
r=zeros(n,n);

prevr=r;    
comparison=50;
maxtrix=zeros(n-1,1);
arrlist=zeros(1,n-1);

while comparison>=1
    for i=1:n
        for k=1:n
            for k1=1:n
                if k1==k
                    continue
                else
                    arrlist(k1)=s(i,k1)+preva(i,k1);
                end
            end
            r(i,k)=(1-lambda)*(s(i,k)-max(arrlist)) + lambda*prevr(i,k);
            end
            end
     for i=1:n
        for k=1:n
            for i1=1:n
                if i1==i || i1==k
                    continue
                else
                    maxtrix(i1)=max(0,r(i1,k));
                end
            end
            smax=sum(maxtrix);
            if k==i
                a(i,k)=(1-lambda)*(smax + max(0,r(k,k))) + lambda* preva(i,k);
            else
                a(i,k)=(1-lambda)*(min(0,r(k,k)+smax))+ lambda * preva(i,k);
            end
            
        end
    end
   comparison=comparison-1;
   preva=a;
   prevr=r; 
    
end
 
exemplar=[];
data=zeros(1,n);

for i=1:n
 Matr=zeros(1,n);
    for j=1:n
      Matr(j)=r(i,j)+a(i,j);
     end
    [p,q]=max(Matr);
     data(i)=q;
     fprintf("i= %d belongs to cluster q=%d\n",i, q);
     if Arr(i,:)==Arr(q,:)
         exemplar(end+1,:)=Arr(q, :);
     end
     
end

%exemplar check
for i=1:size(exemplar,1)
    temp2=find(ismember(Arr,exemplar(i,:)),1);
    data(temp2)=temp2;
end

c=length(exemplar);
U=zeros(c,n);

for i=1:n
     temp=find(exemplar==data(i));
     U(temp,i)=1;
    end

% 
% graph(Arr,exemplar,U);
% 
% 
% function graph(Arr,exemplar,U)
%     color=max((max(U)==U).*linspace(0,1/2,size(U,1))')';
%     hold on
%     scatter(Arr(:,1),Arr(:,2),20, color);
%     scatter(exemplar(:,1),exemplar(:,2),100,'filled');
%     hold off
% end
