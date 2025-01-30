function Y = spacer_unit2(x, x0, charging, np, dx,dsp,  daem, dTT,X_poly, V_T, F, C_St_vol, p_pore, p_IEP, Di_Na, Di_Cl,Di_Na_b, Di_Cl_b, I_max, L_elec, c_di_0,c_con_0,  Alfa, z_Na, z_Cl,E,tau_sp, kk, j_Na_left_flux, j_Na_right_flux,j_Cl_left_flux, j_Cl_right_flux)
    Y = zeros(9*(np+1)+6, 1);
%% concentrate spacer channel 1
    for i = 2:np
        Y(0*(np+1)+i) = Di_Na_b* (   (x(0*(np+1)+i+1)+x(0*(np+1)+i-1)-2*x(0*(np+1)+i))/dsp/dsp + z_Na*(x(0*(np+1)+i+1)-x(0*(np+1)+i-1))/2/dsp*(x(2*(np+1)+i+1)-x(2*(np+1)+i-1))/2/dsp + z_Na*x(0*(np+1)+i)*(x(2*(np+1)+i+1)+x(2*(np+1)+i-1)-2*x(2*(np+1)+i))/dsp/dsp   ) -(x(0*(np+1)+i)-x0(0*(np+1)+i))/dTT + (c_con_0 - x(0*(np+1)+i))/tau_sp;
        Y(1*(np+1)+i) = Di_Cl_b* (   (x(1*(np+1)+i+1)+x(1*(np+1)+i-1)-2*x(1*(np+1)+i))/dsp/dsp + z_Cl*(x(1*(np+1)+i+1)-x(1*(np+1)+i-1))/2/dsp*(x(2*(np+1)+i+1)-x(2*(np+1)+i-1))/2/dsp + z_Cl*x(1*(np+1)+i)*(x(2*(np+1)+i+1)+x(2*(np+1)+i-1)-2*x(2*(np+1)+i))/dsp/dsp   ) -(x(1*(np+1)+i)-x0(1*(np+1)+i))/dTT + (c_con_0 - x(1*(np+1)+i))/tau_sp;
    end

    for i = 1:np+1
        Y(2*(np+1)+i) = x(0*(np+1)+i)-x(1*(np+1)+i);
    end
        Y(0*(np+1)+1)    = -(-Di_Na_b*(3*x(0*(np+1)+1)-4*x(0*(np+1)+2)+x(0*(np+1)+3))/2/dsp - Di_Na_b*z_Na*x(0*(np+1)+1)*(3*x(2*(np+1)+1)-4*x(2*(np+1)+2)+x(2*(np+1)+3))/2/dsp) + x(9*(np+1)+2);% -(-Di_Cl_b*(3*x(1*(np+1)+1)-4*x(1*(np+1)+2)+x(1*(np+1)+3))/2/dsp - Di_Cl_b*z_Cl*x(1*(np+1)+1)*(3*x(2*(np+1)+1)-4*x(2*(np+1)+2)+x(2*(np+1)+3))/2/dsp) - j_Cl_right_cse1;
        Y(0*(np+1)+np+1) = x(2*(np+1)+1);
        Y(1*(np+1)+1)    = -Di_Na_b*(3*x(0*(np+1)+np+1)-4*x(0*(np+1)+np)+x(0*(np+1)+np-1))/2/dsp - Di_Na_b*z_Na*x(0*(np+1)+np+1)*(3*x(2*(np+1)+np+1)-4*x(2*(np+1)+np)+x(2*(np+1)+np-1))/2/dsp - x(9*(np+1)+1);%-Di_Na*(3*x(3*(np+1)+1)-4*x(3*(np+1)+2)+x(3*(np+1)+3))/2/daem - Di_Na*z_Na*x(3*(np+1)+1)*(3*x(5*(np+1)+1)-4*x(5*(np+1)+2)+x(5*(np+1)+3))/2/daem;
        Y(1*(np+1)+np+1) = -Di_Cl_b*(3*x(1*(np+1)+np+1)-4*x(1*(np+1)+np)+x(1*(np+1)+np-1))/2/dsp - Di_Cl_b*z_Cl*x(1*(np+1)+np+1)*(3*x(2*(np+1)+np+1)-4*x(2*(np+1)+np)+x(2*(np+1)+np-1))/2/dsp - x(9*(np+1)+2);%-Di_Cl*(3*x(4*(np+1)+1)-4*x(4*(np+1)+2)+x(4*(np+1)+3))/2/daem - Di_Cl*z_Cl*x(4*(np+1)+1)*(3*x(5*(np+1)+1)-4*x(5*(np+1)+2)+x(5*(np+1)+3))/2/daem ;

%% AEM1
    for i = 2:np
        Y(3*(np+1)+i) = -Di_Na* ( x(3*(np+1)+i+1)-x(3*(np+1)+i-1))/2/daem - z_Na*Di_Na* x(3*(np+1)+i)*(x(5*(np+1)+i+1)-x(5*(np+1)+i-1))/2/daem - x(9*(np+1)+1);
        Y(4*(np+1)+i) =  -Di_Cl* ( x(4*(np+1)+i+1)-x(4*(np+1)+i-1))/2/daem - z_Cl*Di_Cl* x(4*(np+1)+i)*(x(5*(np+1)+i+1)-x(5*(np+1)+i-1))/2/daem - x(9*(np+1)+2);
       
    end

    for i = 1:np+1
        Y(5*(np+1)+i) = x(3*(np+1)+i)-x(4*(np+1)+i)+X_poly;
    end
        Y(3*(np+1)+1)    = x(5*(np+1)+1)-x(2*(np+1)+np+1) + asinh(-X_poly/(2*x(0*(np+1)+np+1)));
        Y(3*(np+1)+np+1) = x(5*(np+1)+np+1)-x(8*(np+1)+1) + asinh(-X_poly/(2*x(6*(np+1)+1)));
        Y(4*(np+1)+1)    = (x(3*(np+1)+1)+x(4*(np+1)+1))^2 - X_poly^2 - (2*x(0*(np+1)+np+1))^2;
        Y(4*(np+1)+np+1) = (x(3*(np+1)+np+1)+x(4*(np+1)+np+1))^2 - X_poly^2 - (2*x(6*(np+1)+1))^2;

        if charging == 1
            Y(9*(np+1)+3)    = x(9*(np+1)+1) - x(9*(np+1)+2) - (j_Na_left_flux - j_Cl_left_flux);
        elseif charging == 2    
            Y(9*(np+1)+3)    = x(9*(np+1)+1) - x(9*(np+1)+2) ;
            
        end
        
        Y(9*(np+1)+1) = - Di_Cl*(3*x(4*(np+1)+np+1)-4*x(4*(np+1)+np)+x(4*(np+1)+np-1))/2/daem - Di_Cl*z_Cl*x(4*(np+1)+np+1)*(3*x(5*(np+1)+np+1)-4*x(5*(np+1)+np)+x(5*(np+1)+np-1))/2/daem - x(9*(np+1)+2);
        Y(9*(np+1)+2) = -(-Di_Cl*(3*x(4*(np+1)+1)-4*x(4*(np+1)+2)+x(4*(np+1)+3))/2/daem - Di_Cl*z_Cl*x(4*(np+1)+1)*(3*x(5*(np+1)+1)-4*x(5*(np+1)+2)+x(5*(np+1)+3))/2/daem) - x(9*(np+1)+2);
%% diluate spacer channel 1
    for i = 2:np
        Y(6*(np+1)+i) = Di_Na_b* (   (x(6*(np+1)+i+1)+x(6*(np+1)+i-1)-2*x(6*(np+1)+i))/dsp/dsp + z_Na*(x(6*(np+1)+i+1)-x(6*(np+1)+i-1))/2/dsp*(x(8*(np+1)+i+1)-x(8*(np+1)+i-1))/2/dsp + z_Na*x(6*(np+1)+i)*(x(8*(np+1)+i+1)+x(8*(np+1)+i-1)-2*x(8*(np+1)+i))/dsp/dsp   ) -(x(6*(np+1)+i)-x0(6*(np+1)+i))/dTT + (c_di_0 - x(6*(np+1)+i))/tau_sp;
        Y(7*(np+1)+i) = Di_Cl_b* (   (x(7*(np+1)+i+1)+x(7*(np+1)+i-1)-2*x(7*(np+1)+i))/dsp/dsp + z_Cl*(x(7*(np+1)+i+1)-x(7*(np+1)+i-1))/2/dsp*(x(8*(np+1)+i+1)-x(8*(np+1)+i-1))/2/dsp + z_Cl*x(7*(np+1)+i)*(x(8*(np+1)+i+1)+x(8*(np+1)+i-1)-2*x(8*(np+1)+i))/dsp/dsp   ) -(x(7*(np+1)+i)-x0(7*(np+1)+i))/dTT + (c_di_0 - x(7*(np+1)+i))/tau_sp;
    end

    for i = 1:np+1
        Y(8*(np+1)+i) = x(6*(np+1)+i)-x(7*(np+1)+i);
    end
        Y(6*(np+1)+1)    = -(-Di_Cl_b*(3*x(7*(np+1)+1)-4*x(7*(np+1)+2)+x(7*(np+1)+3))/2/dsp - Di_Cl_b*z_Cl*x(7*(np+1)+1)*(3*x(8*(np+1)+1)-4*x(8*(np+1)+2)+x(8*(np+1)+3))/2/dsp) - x(9*(np+1)+2);
        Y(6*(np+1)+np+1) = x(9*(np+1)+3) - (-Di_Na_b*(3*x(0*(np+1)+np+1)-4*x(0*(np+1)+np)+x(0*(np+1)+np-1))/2/dsp - Di_Na_b*z_Na*x(0*(np+1)+np+1)*(3*x(2*(np+1)+np+1)-4*x(2*(np+1)+np)+x(2*(np+1)+np-1))/2/dsp);
        Y(7*(np+1)+1)    = -Di_Na_b*(3*x(6*(np+1)+np+1)-4*x(6*(np+1)+np)+x(6*(np+1)+np-1))/2/dsp - Di_Na_b*z_Na*x(6*(np+1)+np+1)*(3*x(8*(np+1)+np+1)-4*x(8*(np+1)+np)+x(8*(np+1)+np-1))/2/dsp + x(9*(np+1)+2);
        Y(7*(np+1)+np+1) = -Di_Cl_b*(3*x(7*(np+1)+np+1)-4*x(7*(np+1)+np)+x(7*(np+1)+np-1))/2/dsp - Di_Cl_b*z_Cl*x(7*(np+1)+np+1)*(3*x(8*(np+1)+np+1)-4*x(8*(np+1)+np)+x(8*(np+1)+np-1))/2/dsp + x(9*(np+1)+1);

        Y(9*(np+1)+4)    = x(9*(np+1)+4) - (-Di_Cl_b*(3*x(1*(np+1)+np+1)-4*x(1*(np+1)+np)+x(1*(np+1)+np-1))/2/dsp - Di_Cl_b*z_Cl*x(1*(np+1)+np+1)*(3*x(2*(np+1)+np+1)-4*x(2*(np+1)+np)+x(2*(np+1)+np-1))/2/dsp);
        Y(9*(np+1)+5)    = x(9*(np+1)+5) - (- (-Di_Na_b*(3*x(6*(np+1)+1)-4*x(6*(np+1)+2)+x(6*(np+1)+3))/2/dsp - Di_Na_b*z_Na*x(6*(np+1)+1)*(3*x(8*(np+1)+1)-4*x(8*(np+1)+2)+x(8*(np+1)+3))/2/dsp));
        Y(9*(np+1)+6)    = x(9*(np+1)+6) - (-(-Di_Cl_b*(3*x(7*(np+1)+1)-4*x(7*(np+1)+2)+x(7*(np+1)+3))/2/dsp - Di_Cl_b*z_Cl*x(7*(np+1)+1)*(3*x(8*(np+1)+1)-4*x(8*(np+1)+2)+x(8*(np+1)+3))/2/dsp));
