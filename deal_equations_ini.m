function Y = deal_equations_ini(x, x0, np, V_T, F, C_St_vol, Alfa, z_Na, z_Cl,E,phi_IEP_ini,c_v_ini, phi_interfacial_ini)
    Y = zeros(4 * (np + 1)+1, 1);
    for i = 1:np+1
        Y(0*(np+1)+i) = x(0*(np+1)+i) - c_v_ini(i)*exp(-z_Na*x(2*(np+1)+i)+E/(x(0*(np+1)+i)+x(1*(np+1)+i)));    
        Y(1*(np+1)+i) = x(1*(np+1)+i) - c_v_ini(i)*exp(-z_Cl*x(2*(np+1)+i)+E/(x(0*(np+1)+i)+x(1*(np+1)+i)));
        Y(2*(np+1)+i) = x(3*(np+1)+i) + (x(0*(np+1)+i)-x(1*(np+1)+i))/(C_St_vol*V_T/F+Alfa*(x(0*(np+1)+i)-x(1*(np+1)+i))^2*V_T/F);
        Y(3*(np+1)+i) = phi_IEP_ini(i) +  phi_interfacial_ini(i) + x(2*(np+1)+i) + x(3*(np+1)+i) - x(4*(np+1)+1);
    end     
    Y(4*(np+1)+1) = z_Na*(x(0*(np+1)+1)+2*sum(x(0*(np+1)+2:0*(np+1)+np))+x(0*(np+1)+np+1))+z_Cl*( x(np+1+1)+2*sum(x(np+1+2:np+1+np))+x(np+1+np+1));

end