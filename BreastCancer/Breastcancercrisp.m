data=importdata('breastCancer.txt');
 x=data(:,1:8);
 out=data(:,9);
 
 n= length(x);
S=zeros(n,n);%initalizing the similarity matrix to zero


%Computing the Similarity Matrix 
for i=1:n
    for k=1:n
        S(i,k)=-norm(x(i)-x(k)); %Similarity matrix is computed by taking negative norm of two points here, however, other methods can be used
    end
end

I=(median(median(S))*eye(n))*170; %Evaluating I using median of the off-diagonal elements of Similarity matrix
S=S+I;





% Initialize messages
N=size(S,1); 
A=zeros(N,N); %Initalizing availability matrix to zero
R=zeros(N,N); %Initalizing responsibility matrix to zero
S=S+(eps*randn(N,N))*(max(S(:))-min(S(:)))*2; % Remove degeneracies
lam=0.92; % Set damping factor


for i=1:100 %number of iterations set to 100, can be changed depending on dataset
  
    % Compute responsibilities
    Rold=R; %Using Rold to record the previous iteration's responsibility matrix
    AS=A+S; 
    [Y,I]=max(AS,[],2); %Y stores the value of maximum value of each column of AS matrix, I stores the row index of the that value in AS
    for i=1:N 
      AS(i,I(i))=-realmax; %setting the value of maximum of AS (where i=k) equal to the IEEE minimum value
    end
    [Y2,I2]=max(AS,[],2); %second maximum value(for i=k) case
    R=S-repmat(Y,[1,N]); %repmat function makes a matrix with all N columns equal to column vector Y. Therefore, R=S-Max(A+S)
    for i=1:N 
      R(i,I(i))=S(i,I(i))-Y2(i); 
    end
    R=(1-lam)*R+lam*Rold; % Dampen responsibilities

    % Compute availabilities
    Aold=A;  %Using Aold to record the previous iteration's availibility matrix
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
idx=I(c); % Assignments, idx shows the indices of exemplars in the input data
disp(idx);
           

% for i=1:n
%     disp([i,idx(i)]);
% end

disp([unique(idx),histc(idx,unique(idx))]); %displays the exemplar index and the number of elements belonging to its cluster


acc=zeros(2,2);
y=zeros(1,n);
ybar=zeros(1,n);

for i=1:n
    if (idx(i)>=1 && idx(i)<=458)
        ybar(i)=2;
    elseif (idx(i)>=459 && idx(i)<=699)
        ybar(i)=4;
    end
end
for i=1:n
    acc(out(i)/2,ybar(i)/2)= acc(out(i)/2,ybar(i)/2)+1;
end

disp(acc);