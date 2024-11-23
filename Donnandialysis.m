function Y = Donnandialysis(x, x0, np, dx, X_poly, Di_Na, Di_Cl, c_di_0,c_con_0,  z_Na, z_Cl)
    Y = zeros(3 * (np + 1)+2, 1);    
    for i = 2:np
        Y(i) = -Di_Na* ( x(0*(np+1)+i+1)-x(0*(np+1)+i-1))/2/dx - z_Na*Di_Na* x(0*(np+1)+i)*(x(2*(np+1)+i+1)-x(2*(np+1)+i-1))/2/dx - x(3*(np+1)+2);
        Y(1*(np+1)+i) =  -Di_Cl* ( x(1*(np+1)+i+1)-x(1*(np+1)+i-1))/2/dx - z_Cl*Di_Cl* x(1*(np+1)+i)*(x(2*(np+1)+i+1)-x(2*(np+1)+i-1))/2/dx - x(3*(np+1)+2);
    end
    for i = 1:np+1
        Y(2*(np+1)+i) = x(0*(np+1)+i)-x(1*(np+1)+i)-X_poly;
    end     
    Y(1) = x(0*(np+1)+1) - c_di_0 *exp(-z_Na*x(2*(np+1)+1));
    Y(np+1) = x(1*(np+1)+1) - c_di_0 *exp(-z_Cl*x(2*(np+1)+1));
    Y(np+2) = c_con_0 - x(0*(np+1)+np+1) *exp(-z_Na*(x(3*(np+1)+1)-x(2*(np+1)+1)));
    Y(2*np+2) = c_con_0 - x(1*(np+1)+np+1) *exp(-z_Cl*(x(3*(np+1)+1)-x(2*(np+1)+1)));  
    Y(3*(np+1)+1) = x(3*(np+1)+2) -Di_Na*z_Na*x(0*(np+1)+1)*(3*x(2*(np+1)+1)-4*x(2*(np+1)+2)+x(2*(np+1)+3))/2/dx-Di_Na*(3*x(0*(np+1)+1)-4*x(0*(np+1)+2)+x(0*(np+1)+3))/2/dx;
    Y(3*(np+1)+2) = x(3*(np+1)+2) +Di_Na*z_Na*x(0*(np+1)+np+1)*(3*x(2*(np+1)+np+1)-4*x(2*(np+1)+np)+x(2*(np+1)+np-1))/2/dx+Di_Na*(3*x(0*(np+1)+np+1)-4*x(0*(np+1)+np)+x(0*(np+1)+np-1))/2/dx;
end
