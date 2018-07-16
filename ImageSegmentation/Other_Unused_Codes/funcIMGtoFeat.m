%function x=getdata(filename1,filename2,filename3,filename4,interval)
n=3;
I=imread('s1.jpg');
%I=rgb2gray(I);
I=double(I);
miny=min(min(I))*ones(size(I));
maxy=max(max(I))*ones(size(I));
img=(I-miny)./(maxy-miny);
map=zeros(size(I));
interval=5;
for i=1:n

    map=map+((i-1)*255/(n-1))*(img(:,:,1)>=((i-1)/(n-1))-(1/(2*(n-1))) & img(:,:,1)<((i-1)/(n-1))+(1/(2*(n-1))));
end
%map=0*(img(:,:,1)>=0 & img(:,:,1)<0.1)+51*(img(:,:,1)>=0.1 & img(:,:,1)<0.3)+102*(img(:,:,1)>=0.3 & img(:,:,1)<0.5)+153*(img(:,:,1)>=0.5 & img(:,:,1)<0.7)+204*(img(:,:,1)>=0.7 & img(:,:,1)<0.9)+255*(img(:,:,1)>=0.9 & img(:,:,1)<=1);
% I=reshape(I,[],1);
check=size(map,1);

indarr=[];
val=zeros(n,1);
 for row=1:interval:floor(size(map,1)-interval)
        for col=1:interval:floor(size(map,2)-interval)
            i=row+floor(rand*interval);
            j=col+floor(rand*interval);
            temp=map(i,j);
            
             for it=1:n
               if temp==(it-1)*255/(n-1) 
                temparr=(j-1)*check+i;
                indarr(end+1)=temparr;
                val(it,1)=val(it,1)+1;
               end 
            end
%             if temp==0 && val(1,1)<100
%                 temparr=(j-1)*check+i;
%                 indarr(end+1)=temparr;
%                 val(1,1)=val(1,1)+1;
%             end
%             if temp==51 && val(2,1)<100
%                temparr=(j-1)*check+i;
%                 indarr(end+1)=temparr;
%                 val(2,1)=val(2,1)+1;
%             end
%             if temp==102 && val(3,1)<100
%                 temparr=(j-1)*check+i;
%                 indarr(end+1)=temparr;
%                 val(3,1)=val(3,1)+1;
%             end
%             if temp==153 && val(4,1)<100
%                 temparr=(j-1)*check+i;
%                 indarr(end+1)=temparr;
%                 val(4,1)=val(4,1)+1;
%             end
%             if temp==204 && val(5,1)<100
%                 temparr=(j-1)*check+i;
%                 indarr(end+1)=temparr;
%                 val(5,1)=val(5,1)+1;
%             end
%             if temp==255 && val(6,1)<100
%                 temparr=(j-1)*check+i;
%                 indarr(end+1)=temparr;
%                 val(6,1)=val(6,1)+1;
%             end
            
        end
 end
 disp(val)
 
 I1=getfeatures('s1_entro.jpg',indarr);
 I2=getfeatures('s1_med.jpg',indarr);
 I3=getfeatures('s1_gauss.jpg',indarr);
 Out=map(indarr(1,:));
  Out=Out';
 
 x=[I1 I2 I3 Out];
 [dummy index]=sort(x(:,4));
 x=x(index,:);
 



 function I=getfeatures(filename,arr)
  i1=imread(filename);
  i1=double(i1);
  I=i1(arr(1,:));
  I=I';
 end