function Y = deal_equations(x, x0, np, dx, dTT,X_poly, V_T, F, C_St_vol, p_pore, p_IEP, Di_Na, Di_Cl, I_max, L_elec, charging_discharge, c_di_0,c_con_0,  Alfa, z_Na, z_Cl,E)
    Y = zeros(9 * (np + 1)+7, 1);
    for i = 2:np
        Y(i) = Di_Na* (   (x(0*(np+1)+i+1)+x(0*(np+1)+i-1)-2*x(0*(np+1)+i))/dx/dx + z_Na*(x(0*(np+1)+i+1)-x(0*(np+1)+i-1))/2/dx*(x(5*(np+1)+i+1)-x(5*(np+1)+i-1))/2/dx + z_Na*x(0*(np+1)+i)*(x(5*(np+1)+i+1)+x(5*(np+1)+i-1)-2*x(5*(np+1)+i))/dx/dx   ) -p_IEP*(x(0*(np+1)+i)-x0(0*(np+1)+i))/dTT-p_pore*(x(2*(np+1)+i)-x0(2*(np+1)+i))/dTT;
        Y(1*(np+1)+i) = Di_Cl* (   (x(1*(np+1)+i+1)+x(1*(np+1)+i-1)-2*x(1*(np+1)+i))/dx/dx + z_Cl*(x(1*(np+1)+i+1)-x(1*(np+1)+i-1))/2/dx*(x(5*(np+1)+i+1)-x(5*(np+1)+i-1))/2/dx + z_Cl*x(1*(np+1)+i)*(x(5*(np+1)+i+1)+x(5*(np+1)+i-1)-2*x(5*(np+1)+i))/dx/dx   ) -p_IEP*(x(1*(np+1)+i)-x0(1*(np+1)+i))/dTT-p_pore*(x(3*(np+1)+i)-x0(3*(np+1)+i))/dTT;
    end

    for i = 1:np+1
        Y(2*(np+1)+i) = x(0*(np+1)+i)-x(1*(np+1)+i)-X_poly;
        Y(3*(np+1)+i) = x(4*(np+1)+i) - x(0*(np+1)+i)*exp(-z_Na*x(6*(np+1)+i));
        Y(4*(np+1)+i) = x(2*(np+1)+i) - x(4*(np+1)+i)*exp(-z_Na*x(7*(np+1)+i)+E/(x(2*(np+1)+i)+x(3*(np+1)+i)));
        Y(5*(np+1)+i) = x(4*(np+1)+i) - x(1*(np+1)+i)*exp(-z_Cl*x(6*(np+1)+i));
        Y(6*(np+1)+i) = x(3*(np+1)+i) - x(4*(np+1)+i)*exp(-z_Cl*x(7*(np+1)+i)+E/(x(2*(np+1)+i)+x(3*(np+1)+i)));
        Y(7*(np+1)+i) = x(8*(np+1)+i) + (x(2*(np+1)+i)-x(3*(np+1)+i))/(C_St_vol*V_T/F+Alfa*(x(2*(np+1)+i)-x(3*(np+1)+i))^2*V_T/F);
        Y(8*(np+1)+i) = x(5*(np+1)+i) + x(6*(np+1)+i) + x(7*(np+1)+i) + x(8*(np+1)+i)-x(9*(np+1)+1);
    end

    if charging_discharge == 1
        Y(0*(np+1)+1) = x(5*(np+1)+1) + asinh(X_poly/(2*c_di_0));
        Y(0*(np+1)+np+1) = x(4*(np+1)+1) - c_di_0;
        Y(1*(np+1)+1) = -F*(Di_Na*x(0*(np+1)+np+1)+Di_Cl*x(1*(np+1)+np+1))*(3*x(5*(np+1)+np+1)-4*x(5*(np+1)+np)+x(5*(np+1)+np-1))-F*(Di_Na-Di_Cl)*(3*x(0*(np+1)+np+1)-4*x(0*(np+1)+np)+x(0*(np+1)+np-1)); 
        Y(1*(np+1)+np+1) = x(4*(np+1)+np+1) - c_con_0;
        Y(9*(np+1)+1) = z_Na*(x(2*(np+1)+1)+2*sum(x(2*(np+1)+2:2*(np+1)+np))+x(2*(np+1)+np+1))+z_Cl*( x(2*(np+1)+np+1+1)+2*sum(x(2*(np+1)+np+1+2:2*(np+1)+np+1+np))+x(2*(np+1)+np+1+np+1)) - x(9*(np+1)+3)*2*np;
        Y(9*(np+1)+2) =  x(9*(np+1)+3)-x0(9*(np+1)+3)-I_max/F*dTT/L_elec/p_pore;  
        Y(9*(np+1)+3) =  p_IEP*(x(1)+2*sum(x(2:np))+x(np+1)+ x(np+1+1)+2*sum(x(np+1+2:np+1+np))+x(np+1+np+1)) + p_pore*(x(2*(np+1)+1)+2*sum(x(2*(np+1)+2:2*(np+1)+np))+x(2*(np+1)+np+1)+ x(2*(np+1)+np+1+1)+2*sum(x(2*(np+1)+np+1+2:2*(np+1)+np+1+np))+x(2*(np+1)+np+1+np+1)) - (p_IEP*(x0(1)+2*sum(x0(2:np))+x0(np+1)+ x0(np+1+1)+2*sum(x0(np+1+2:np+1+np))+x0(np+1+np+1)) + p_pore*(x0(2*(np+1)+1)+2*sum(x0(2*(np+1)+2:2*(np+1)+np))+x0(2*(np+1)+np+1)+ x0(2*(np+1)+np+1+1)+2*sum(x0(2*(np+1)+np+1+2:2*(np+1)+np+1+np))+x0(2*(np+1)+np+1+np+1))) - x(9*(np+1)+2)*dTT/L_elec*2*np;
        Y(9*(np+1)+4) = x(9*(np+1)+4) -Di_Na*z_Na*x(0*(np+1)+1)*(3*x(5*(np+1)+1)-4*x(5*(np+1)+2)+x(5*(np+1)+3))/2/dx-Di_Na*(3*x(0*(np+1)+1)-4*x(0*(np+1)+2)+x(0*(np+1)+3))/2/dx;
        Y(9*(np+1)+5) = x(9*(np+1)+5) +Di_Na*z_Na*x(0*(np+1)+np+1)*(3*x(5*(np+1)+np+1)-4*x(5*(np+1)+np)+x(5*(np+1)+np-1))/2/dx+Di_Na*(3*x(0*(np+1)+np+1)-4*x(0*(np+1)+np)+x(0*(np+1)+np-1))/2/dx;
        Y(9*(np+1)+6) = x(9*(np+1)+6) -Di_Cl*z_Cl*x(1*(np+1)+1)*(3*x(5*(np+1)+1)-4*x(5*(np+1)+2)+x(5*(np+1)+3))/2/dx-Di_Cl*(3*x(1*(np+1)+1)-4*x(1*(np+1)+2)+x(1*(np+1)+3))/2/dx;
        Y(9*(np+1)+7) = x(9*(np+1)+7) +Di_Cl*z_Cl*x(1*(np+1)+np+1)*(3*x(5*(np+1)+np+1)-4*x(5*(np+1)+np)+x(5*(np+1)+np-1))/2/dx+Di_Cl*(3*x(1*(np+1)+np+1)-4*x(1*(np+1)+np)+x(1*(np+1)+np-1))/2/dx;
     elseif charging_discharge ==2
       
        Y(0*(np+1)+1) = x(5*(np+1)+np+1) + asinh(X_poly/(2*c_con_0));
        Y(0*(np+1)+np+1) = x(4*(np+1)+1) - c_di_0;
        Y(1*(np+1)+1) = -F*(Di_Na*x(0*(np+1)+1)+Di_Cl*x(1*(np+1)+1))*(3*x(5*(np+1)+1)-4*x(5*(np+1)+2)+x(5*(np+1)+3))-F*(Di_Na-Di_Cl)*(3*x(0*(np+1)+1)-4*x(0*(np+1)+2)+x(0*(np+1)+3));  
        Y(1*(np+1)+np+1) = x(4*(np+1)+np+1) - c_con_0;

        Y(9*(np+1)+1) = z_Na*(x(2*(np+1)+1)+2*sum(x(2*(np+1)+2:2*(np+1)+np))+x(2*(np+1)+np+1))+z_Cl*( x(2*(np+1)+np+1+1)+2*sum(x(2*(np+1)+np+1+2:2*(np+1)+np+1+np))+x(2*(np+1)+np+1+np+1)) - x(9*(np+1)+3)*2*np;
        Y(9*(np+1)+2) =  x(9*(np+1)+3)-x0(9*(np+1)+3)+I_max/F*dTT/L_elec/p_pore;  
        Y(9*(np+1)+3) =  p_IEP*(x(1)+2*sum(x(2:np))+x(np+1)+ x(np+1+1)+2*sum(x(np+1+2:np+1+np))+x(np+1+np+1)) + p_pore*(x(2*(np+1)+1)+2*sum(x(2*(np+1)+2:2*(np+1)+np))+x(2*(np+1)+np+1)+ x(2*(np+1)+np+1+1)+2*sum(x(2*(np+1)+np+1+2:2*(np+1)+np+1+np))+x(2*(np+1)+np+1+np+1)) - (p_IEP*(x0(1)+2*sum(x0(2:np))+x0(np+1)+ x0(np+1+1)+2*sum(x0(np+1+2:np+1+np))+x0(np+1+np+1)) + p_pore*(x0(2*(np+1)+1)+2*sum(x0(2*(np+1)+2:2*(np+1)+np))+x0(2*(np+1)+np+1)+ x0(2*(np+1)+np+1+1)+2*sum(x0(2*(np+1)+np+1+2:2*(np+1)+np+1+np))+x0(2*(np+1)+np+1+np+1))) - x(9*(np+1)+2)*dTT/L_elec*2*np;
        Y(9*(np+1)+4) = x(9*(np+1)+4) -Di_Na*z_Na*x(0*(np+1)+1)*(3*x(5*(np+1)+1)-4*x(5*(np+1)+2)+x(5*(np+1)+3))/2/dx-Di_Na*(3*x(0*(np+1)+1)-4*x(0*(np+1)+2)+x(0*(np+1)+3))/2/dx;
        Y(9*(np+1)+5) = x(9*(np+1)+5) +Di_Na*z_Na*x(0*(np+1)+np+1)*(3*x(5*(np+1)+np+1)-4*x(5*(np+1)+np)+x(5*(np+1)+np-1))/2/dx+Di_Na*(3*x(0*(np+1)+np+1)-4*x(0*(np+1)+np)+x(0*(np+1)+np-1))/2/dx;
        Y(9*(np+1)+6) = x(9*(np+1)+6) -Di_Cl*z_Cl*x(1*(np+1)+1)*(3*x(5*(np+1)+1)-4*x(5*(np+1)+2)+x(5*(np+1)+3))/2/dx-Di_Cl*(3*x(1*(np+1)+1)-4*x(1*(np+1)+2)+x(1*(np+1)+3))/2/dx;
        Y(9*(np+1)+7) = x(9*(np+1)+7) +Di_Cl*z_Cl*x(1*(np+1)+np+1)*(3*x(5*(np+1)+np+1)-4*x(5*(np+1)+np)+x(5*(np+1)+np-1))/2/dx+Di_Cl*(3*x(1*(np+1)+np+1)-4*x(1*(np+1)+np)+x(1*(np+1)+np-1))/2/dx;
    end
end
