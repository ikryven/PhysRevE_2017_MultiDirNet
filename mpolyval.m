function z=mpolyval(U,x,y)

p=zeros(length(x),size(U,2));
for i=1:size(U,2)
    
    p(1:length(x),i) = polyval( U(:,i), x );
    
end;

for i=1:length(x)
    z(i) = polyval(p(i,:),y(i));
end;