function [reg1,reg2, MSD,minimum]=bestreg(y,q1,q2,p,H,w_true,lambda1, lambda2) 

        num1=length(lambda1);        
        num2=length(lambda2);
        MSD=zeros(num1,num2);        
  for i=1:num1
      for j=1:num2
                MSD(i,j)=JSLMSP2(y,q1,q2,p,H,w_true,lambda1(i),lambda2(j));
      end
  end
    minimum=min(min(MSD));
    [row,col]=find(MSD==minimum);
    reg1=lambda1(row);
    reg2=lambda2(col);    
end