function b=array_form(element_array,d,lam,theta,alfa,tag)   % element_arrayΪ��Դ�� dΪ��Դ��� lamΪ������thetaΪ����ǣ�alfaΪ��λ�� tagΪ��Դģ��
%%%%%%% ��Ԫ��������%%%%%%
switch tag
    case 1 %������
        mm=0:1:(element_array-1);
        a=d;
        a=sqrt(2)/2*d;
        am=exp(-1j*2*pi*a/lam*sin(theta*pi/180)*cos(alfa*pi/180*ones(1,element_array)-(2*mm+1)*pi/element_array));
        
     case 2 %����
        mm=0:1:(element_array-1);
        a=sqrt(2)/2*d;
        am=exp(-1j*2*pi*a/lam*sin(theta*pi/180)*cos(alfa*pi/180*ones(1,element_array)-2*mm*pi/element_array));
        
     case 3 %y��
        mm=0:1:(element_array-2);
        a=d;
        am=exp(-1j*2*pi*a/lam*sin(theta*pi/180)*cos(alfa*pi/180*ones(1,element_array-1)-(3*mm+1)*pi/element_array));
        am=[1 am];
        
     case 4 %Բ�󣬲�����������
        mm=0:1:(element_array-1);
        a=d;
        am=exp(-1j*2*pi*a/lam*sin(theta*pi/180)*cos(alfa*pi/180*ones(1,element_array)-2*mm*pi/element_array));
      
     case 5 %Բ�󣬺���������
        mm=0:1:(element_array-2);
        a = d;
        am=exp(-1j*2*pi*a/lam*sin(theta*pi/180)*cos(alfa*pi/180*ones(1,element_array-1)-2*mm*pi/element_array));
        am=[1 am] ;
        
      case 6 %����
        mm=0:1:(element_array-1);
        a = d;
        am=exp(-1j*2*pi*a/lam*sin(theta*pi/180)*mm');
        %am=[1 am] ;  
end
b=am.';