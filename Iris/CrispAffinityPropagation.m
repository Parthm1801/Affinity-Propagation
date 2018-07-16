%n is the number of data point
%epsilon is the threshold of stopping algorithm with regard to membership value of the datapoints
%lam is the damping factor
%S is a n*n simlarity matrix, where n is the number of datapoints
%I is a n*n diagonal matrix, used to evaluate the diagonal elements of the Similarity Matrix, also called preferences
%A is a n*n availability matrix
%R is a n*n responsibility matrix

%Input/Dataset
%x=[1 2 3 4 5 6 7  997 998 996  999 4000 4002 4004];
%x= [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 52 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 79 80 81 82 83 84 85 86 88 89 90 91 92 93 94 95 96 97 98 99 100 102 107 114 120 122 127 128 139 143  51 53 77 78 87 101 103 104 105 106 108 109 110 111 112 113 115 116 117 118 119 121 123 124 125 126 129 130 131 132 133 134 135 136 137 138 140 141 142 144 145 146 147 148 149 150];
%x=[3 4 3 2 1; 4 3 5 1 1;3 5 3 3 3; 2 1 3 3 2;1 1 3 2 3];
%x=[1 2 4 5; 4 5 7 8 ;4 6 6 7; 99 99 99 99; 99 99 98 91];
%x=[1 2 3 91 99 4000 4005]';
% delimiterIn =',';
data=split(importdata('iris dataset.txt'),",");
size(data)
% disp(data)
X=data(:,1:4);
out=data(:,5);
% disp(X)
% disp(Y)
x=str2double(X);
%x=sort(x);


n= length(x);
S=zeros(n,n);%initalizing the similarity matrix to zero


%Computing the Similarity Matrix 
for i=1:n
    for k=1:n
        S(i,k)=-norm(x(i)-x(k)); %Similarity matrix is computed by taking negative norm of two points here, however, other methods can be used
    end
end

I=(median(median(S))*eye(n))/0.05; %Evaluating I using median of the off-diagonal elements of Similarity matrix
S=S+I;





% Initialize messages
N=size(S,1); 
A=zeros(N,N); %Initalizing availability matrix to zero
R=zeros(N,N); %Initalizing responsibility matrix to zero
S=S+(eps*randn(N,N))*(max(S(:))-min(S(:))); % Remove degeneracies
lam=0.9; % Set damping factor


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

            

for i=1:n
    disp([i,idx(i)]);
end

disp([unique(idx),histc(idx,unique(idx))]); %displays the exemplar index and the number of elements belonging to its cluster


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

color=(idx)/(sum(uidx));
%color=[color color color]+[rand() rand() rand()]
scatter3((x(:,2)-min(x(:,2)))/(max(x(:,2))-min(x(:,2))),(x(:,3)-min(x(:,3)))/(max(x(:,3))-min(x(:,3))),(x(:,4)-min(x(:,4)))/(max(x(:,4))-min(x(:,4))),100,color,'filled');