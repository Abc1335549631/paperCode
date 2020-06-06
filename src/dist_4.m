function [Sum]=dist_4(x)

y=x;
D =[];
n=size(x,1);
for  i=1:n
    for j=1:n
     d=norm(x(i,:)-y(j,:));
     D =[D d];
        
    end     
end       
  Sum=(sum(D ) / ( n*n  ));

end

