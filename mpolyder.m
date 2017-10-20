function dU=mpolyder(U,d)



dU = zeros(size(U));

if d==1
    
for i=1:size(U,2)
    der=polyder(U(:,i));
    dU(end-length(der)+1:end,i) =der;
    
end;
end;


if d==2

for i=1:size(U,1)
    
    der=polyder(U(i,:));
    
    dU(i,end-length(der)+1:end) = der;
    
end;

end;

if all( dU(:,1) ==0 )
    dU(:,1)=[];
end

if all( dU(1,:) ==0 )
    dU(1,:)=[];
end