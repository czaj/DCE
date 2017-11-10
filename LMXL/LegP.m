function L = LegP(Order, err)

[A,B] = size(err);

L = zeros(A,B, Order);
L(:,:,1) = err;
 L(:,:,2) = ((2*2-1)/2)*err.*L(:,:,1)-(2-1)/2;
 L(:,:,3) = ((2*3-1)/3)*err.*L(:,:,2)-((3-1)/3)*L(:,:,1);
 if Order > 3
     for i =1:(Order-3)
          Tmp = 3+i;
          L(:,:, Tmp) = ((2*Tmp-1)/Tmp)*err.*L(:,:,Tmp-1)-((Tmp-1)/Tmp)*L(:,:,Tmp-2);
     end
 end
