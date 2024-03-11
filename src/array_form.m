function b=array_form(element_array,d,lam,theta,alfa,tag)   % element_array为阵源数 d为阵源间距 lam为波长，theta为入射角，alfa为方位角 tag为阵源模型
%%%%%%% 四元阵列流型%%%%%%
switch tag
    case 1 %正方形
        mm=0:1:(element_array-1);
        a=d;
        a=sqrt(2)/2*d;
        am=exp(-1j*2*pi*a/lam*sin(theta*pi/180)*cos(alfa*pi/180*ones(1,element_array)-(2*mm+1)*pi/element_array));
        
     case 2 %棱形
        mm=0:1:(element_array-1);
        a=sqrt(2)/2*d;
        am=exp(-1j*2*pi*a/lam*sin(theta*pi/180)*cos(alfa*pi/180*ones(1,element_array)-2*mm*pi/element_array));
        
     case 3 %y形
        mm=0:1:(element_array-2);
        a=d;
        am=exp(-1j*2*pi*a/lam*sin(theta*pi/180)*cos(alfa*pi/180*ones(1,element_array-1)-(3*mm+1)*pi/element_array));
        am=[1 am];
        
     case 4 %圆阵，不含中心阵子
        mm=0:1:(element_array-1);
        a=d;
        am=exp(-1j*2*pi*a/lam*sin(theta*pi/180)*cos(alfa*pi/180*ones(1,element_array)-2*mm*pi/element_array));
      
     case 5 %圆阵，含中心阵子
        mm=0:1:(element_array-2);
        a = d;
        am=exp(-1j*2*pi*a/lam*sin(theta*pi/180)*cos(alfa*pi/180*ones(1,element_array-1)-2*mm*pi/element_array));
        am=[1 am] ;
        
      case 6 %线阵
        mm=0:1:(element_array-1);
        a = d;
        am=exp(-1j*2*pi*a/lam*sin(theta*pi/180)*mm');
        %am=[1 am] ;  
end
b=am.';