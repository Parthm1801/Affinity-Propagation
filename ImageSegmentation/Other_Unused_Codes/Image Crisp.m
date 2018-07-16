x=get_image_data('eagle.jpg',5);

n= size(x,1);
S=zeros(n,n);%initalizing the similarity matrix to zero


%Computing the Similarity Matrix 
for i=1:n
    for k=1:n
        S(i,k)=-norm(x(i,:)-x(k,:)); %Similarity matrix is computed by taking negative norm of two points here, however, other methods can be used
    end
end

I=(median(median(S))*eye(n))*10; %Evaluating I using median of the off-diagonal elements of Similarity matrix
S=S+I;





% Initialize messages
N=size(S,1); 
A=zeros(N,N); %Initalizing availability matrix to zero
R=zeros(N,N); %Initalizing responsibility matrix to zero
S=S+(eps*randn(N,N))*(max(S(:))-min(S(:)))*2; % Remove degeneracies
lam=0.8; % Set damping factor


for iter=1:50 %number of iterations set to 100, can be changed depending on dataset
  
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


% acc=zeros(2,2);
% y=zeros(1,n);
% ybar=zeros(1,n);
% 
% for i=1:n
%     if (idx(i)>=1 && idx(i)<=102)
%         ybar(i)=3;
%     elseif (idx(i)>=103 && idx(i)<=232)
%         ybar(i)=6;
%     end
% end
% for i=1:n
%     acc(out(i)/3,ybar(i)/3)= acc(out(i)/3,ybar(i)/3)+1;
% end
% 
% disp(acc);

uidx=unique(idx);
show_image(x,idx,uidx,'eagle.jpg',5)

    
    
    function X=get_image_data(filename,interval)
    I=imread(filename);
    I_dash=I;
    X=[];
    for row=1:interval:floor(size(I,1)-interval)
        for column=1:interval:floor(size(I,2)-interval)
            x_unfolded=I(row:row+interval-1,column:column+interval-1,:);
            X=[X ;x_unfolded(:)'];
        end
    end
    X=double(X);
end
function show_image(X,idx,uidx,filename,interval)
    color=idx/sum(uidx)  ; 
    color=color*256;
    unique_color=unique(color);
    store=zeros(size(color,1),size(unique_color,1));
    for unique_color_index=1:size(unique_color,1)
        store(:,unique_color_index)=(unique_color(unique_color_index)==color);
    end
   [temp sorted_indices]=sort(sum(store));
   store=store(:,sorted_indices);
   store=store.*unique_color';
   color=sum(store')';
    for i=1:size(X,1)
        X(i,:)=color(i);
    end
    I=imread(filename);
   
    
    data_number=0;
    for row=1:interval:floor(size(I,1)-interval)
        for column=1:interval:floor(size(I,2)-interval)
            data_number=data_number+1;
            x_folded=X(data_number,:);
            I(row:row+interval-1,column:column+interval-1,:)=reshape(x_folded,interval,interval,3);
        end
    end
    
    disp("show");
    imshow(I);
end
   