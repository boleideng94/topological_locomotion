function dy = ODE_func(T,y)

global Nx Ny ks kt kb kl lrb mass J
global damp_theta damp_u Time thst
global actuator_tube progress
global friction_dt friction_time friction_tube 
global mass_act cycle_number_current cycle_number
global y_0 seg lp la
seg = 0; step = 0.01;
if T/Time > progress + 0.0
    progress = progress + step;
    disp([num2str((progress-step)/step),'%']);
    seg = 1;
end
% 
if cycle_number_current == cycle_number
    cycle_number_current = cycle_number_current + 1;
    y_0 = y;
end

if T>friction_time
    flag_fric = 1;
    friction_time = friction_time + friction_dt;
else
    flag_fric = 0;
end

flag_U = 0; flag_V = 0; flag_Th = 0;
% Initial and boundary conditions
dy = zeros(Nx*Ny*8,1);
F = zeros(Nx*Ny*8,1);

% Main part of ODE
for i=1:Ny
    for j=1:Nx
        % i-th row and j-th column
        % Each cell have 3 DOFs: U V Th
        % Define coodinates
        posi = ((i-1)*Nx + j-1)*8;
        up_posi = ((i-1+1)*Nx + j-1)*8;
        down_posi = ((i-1-1)*Nx + j-1)*8;
        left_posi = ((i-1)*Nx + j-1-1)*8;
        right_posi = ((i-1)*Nx + j-1+1)*8;
        %Coordinates of the center mass
        U_cen=posi+1;
        V_cen=posi+2;
        Th_cen=posi+3;
        Rot_cen=posi+4;
        %Coordinates of the top mass
        U_up=up_posi+1;
        V_up=up_posi+2;
        Th_up=up_posi+3;
        Rot_up=up_posi+4;
        %Coordinates of the bottom mass
        U_down=down_posi+1;
        V_down=down_posi+2;
        Th_down=down_posi+3;
        Rot_down=down_posi+4;
        %Coordinates of the right mass
        U_right=right_posi+1;
        V_right=right_posi+2;
        Th_right=right_posi+3;
        Rot_right=right_posi+4;
        %Coordinates of the left mass
        U_left=left_posi+1;
        V_left=left_posi+2;
        Th_left=left_posi+3;
        Rot_left=left_posi+4;
        % force
        force_up = [0,0]; force_down = [0,0]; force_left = [0,0]; force_right = [0,0];
        moment_up = 0; moment_down = 0; moment_left = 0; moment_right = 0;
        shear_on = 1; sig = 1;
        rot_up = 0; rot_down = 0; rot_right = 0; rot_left = 0;

        if j~=Nx % otherwise nothing on the right       
            dU = [y(U_right)-y(U_cen);y(V_right)-y(V_cen)]+((vector(y(Th_right),y(Rot_right),'left')-vector(thst,0,'left')))...
                -((vector(y(Th_cen),y(Rot_cen),'right')-vector(thst,0,'right')));
            dUx = dU(1); dUy = dU(2);
            norm_strain = dUx*cos(y(Rot_cen)) + dUy*sin(y(Rot_cen));
            shear_strain = -dUx*sin(y(Rot_cen)) + dUy*cos(y(Rot_cen));
            
            force_right(1) = force_right(1) + 2*kl*(norm_strain);
            force_right(2) = force_right(2) + 2*ks*(shear_strain);
%             resi_norm = norm_strain-(y(U_right)-y(U_cen) - 4*la/2*(cos(y(Th_right))+cos(y(Th_cen))-2*cos(thst)))
%             resi_shear = shear_strain - (y(V_right)-y(V_cen))
            moment_right = ks*4*la/2*cos(y(Th_cen))*4*la/4*(sin(y(Th_right))-sin(y(Th_cen))) ...
                -force_right(1)*2*4*la/4*sin(y(Th_cen));
            rot_right = y(Rot_right)-y(Rot_cen);
        end
        if i~=Ny %otherwise nothing on the top
            dU = [y(U_up)-y(U_cen);y(V_up)-y(V_cen)]+((vector(y(Th_up),y(Rot_up),'down')-vector(thst,0,'down')))...
                -((vector(y(Th_cen),y(Rot_cen),'up')-vector(thst,0,'up')));
            dUx = dU(1); dUy = dU(2);
%             dU -  [y(U_up)-y(U_cen);y(V_up)-y(V_cen)]
            norm_strain = -(dUx)*sin(y(Rot_cen))+(dUy)*cos(y(Rot_cen));
            shear_strain = (dUx)*cos(y(Rot_cen))+(dUy)*sin(y(Rot_cen));  
            force_up(1) = force_up(1) + 2*ks*(shear_strain);
            force_up(2) = force_up(2) + 2*kl*(norm_strain);
%             resi_norm = norm_strain-(y(V_up)-y(V_cen)-4*la/2*(cos(y(Th_up))+cos(y(Th_cen))-2*cos(thst)))
%             resi_shear = shear_strain - (y(U_up)-y(U_cen))
            moment_up = ks*4*la/2*cos(y(Th_cen))*4*la/4*(sin(y(Th_up))-sin(y(Th_cen))) ...
                -force_up(2)*2*4*la/4*sin(y(Th_cen));
%             rot_up = atan(-(y(U_up)-y(U_cen))/(y(V_up)-y(V_cen)+a));
            rot_up = y(Rot_up)-y(Rot_cen);
        end
        if j~=1 %otherwise nothing on the left
            dU = [y(U_left)-y(U_cen);y(V_left)-y(V_cen)]+((vector(y(Th_left),y(Rot_left),'right')-vector(thst,0,'right')))...
                -((vector(y(Th_cen),y(Rot_cen),'left')-vector(thst,0,'left')));
%             dU
            dUx = dU(1); dUy = dU(2);
            norm_strain = dUx*cos(y(Rot_cen)) + dUy*sin(y(Rot_cen));
            shear_strain = -dUx*sin(y(Rot_cen)) + dUy*cos(y(Rot_cen));
            force_left(1) = force_left(1) + 2*kl*(norm_strain);
            force_left(2) = force_left(2) + 2*ks*(shear_strain);
%             resi_norm = norm_strain-(y(U_left)-y(U_cen)+4*la/2*(cos(y(Th_left))+cos(y(Th_cen))-2*cos(thst)))
%             resi_shear = shear_strain - (y(V_left)-y(V_cen))
            moment_left = ks*4*la/2*cos(y(Th_cen))*4*la/4*(sin(y(Th_left))-sin(y(Th_cen))) ...
                +force_left(1)*2*4*la/4*sin(y(Th_cen));
%             rot_left = atan(-(y(V_left)-y(V_cen))/(-y(U_left)+y(U_cen)+a));
            rot_left = y(Rot_left)-y(Rot_cen);
        end
        if i~=1 %otherwise nothing on the bottom
            dU = [y(U_down)-y(U_cen);y(V_down)-y(V_cen)]+((vector(y(Th_down),y(Rot_down),'up')-vector(thst,0,'up')))...
                -((vector(y(Th_cen),y(Rot_cen),'down')-vector(thst,0,'down')));
            dUx = dU(1); dUy = dU(2);
%             dU - [y(U_down)-y(U_cen);y(V_down)-y(V_cen)]
            norm_strain = -(dUx)*sin(y(Rot_cen))+(dUy)*cos(y(Rot_cen));
            shear_strain = (dUx)*cos(y(Rot_cen))+(dUy)*sin(y(Rot_cen));  
            force_down(1) = force_down(1) + 2*ks*(shear_strain);
            force_down(2) = force_down(2) + 2*kl*(norm_strain);           
%             resi_norm = norm_strain - (y(V_down)-y(V_cen)+4*la/2*(cos(y(Th_down))+cos(y(Th_cen))-2*cos(thst)))
%             resi_shear = shear_strain - (y(U_down)-y(U_cen))
            
            moment_down = ks*4*la/2*cos(y(Th_cen))*4*la/4*(sin(y(Th_down))-sin(y(Th_cen))) ...
                +force_down(2)*2*4*la/4*sin(y(Th_cen));
%             rot_down = atan((y(U_down)-y(U_cen))/(-y(V_down)+y(V_cen)+a));
            rot_down = y(Rot_down)-y(Rot_cen);
        end
       
%         force_up = [0,0]; force_down = [0,0]; 
%         force_left = [0,0]; force_right = [0,0];
        local_moment = 4*kb*lp*(lp*cos(y(Th_cen))-lrb)*sin(y(Th_cen)) - 16*kt*y(Th_cen);
        
        % actuation
        if sum(i == actuator_tube(:,1) & j == actuator_tube(:,2))
            % the element is actuated
            [~,M_input] = func_actuation(y,T,i,j);
        else
            % the element is not actuated
            M_input = 0;
        end
        
        if flag_fric == 1 % update friction
            posi = ((i-1)*Nx + j-1)*8;
            y_fric = y(posi+1:posi+8);
            friction_tube(i,j,:) = func_friction(y_fric,T,i,j);
        end
        % Total force in U & V directions

        F(U_cen) = (force_right(1) + force_left(1) + force_up(1) + force_down(1))*cos(y(Rot_cen))...
                 - (force_right(2) + force_left(2) + force_up(2) + force_down(2))*sin(y(Rot_cen)) ...
                   + friction_tube(i,j,1);
        F(V_cen) = (force_right(2) + force_left(2) + force_up(2) + force_down(2))*cos(y(Rot_cen)) ...
                 + (force_right(1) + force_left(1) + force_up(1) + force_down(1))*sin(y(Rot_cen)) ...
                 + friction_tube(i,j,2);
             
        % limit maximum forces
%         if abs(F(U_cen)>10)
%             'check U'
%             F(U_cen) = 10*sign(F(U_cen));
%         end
%         if abs(F(V_cen)>10)
%             'check V'
%             F(V_cen) = 10*sign(F(V_cen));
%         end
             
        F(Th_cen) = moment_right+moment_left+moment_up+moment_down+local_moment+M_input + 1*friction_tube(i,j,3);
%         F(Rot_cen) = ((rot_up+rot_down+rot_right+rot_left)/4-y(Rot_cen)).*a/2*ks;
        F(Rot_cen) = (force_right(2)-force_up(1)-force_left(2)+force_down(1)).*4*la/2+...
                    (rot_right+rot_up+rot_left+rot_down).*4*la/2*kl;
        
        dy(U_cen) = y(U_cen+4);
        dy(V_cen) = y(V_cen+4);
        dy(Th_cen) = y(Th_cen+4);
        dy(Rot_cen) = y(Rot_cen+4);
        
        if sum(i == actuator_tube(:,1) & j == actuator_tube(:,2))
            mass_total = mass + mass_act;
        else
            mass_total = mass;
        end
        
        dy(U_cen+4) = F(U_cen)/mass_total - damp_u/mass_total*y(U_cen+4);
        dy(V_cen+4) = F(V_cen)/mass_total - damp_u/mass_total*y(V_cen+4);
        dy(Th_cen+4) = F(Th_cen)/J - damp_theta/J*y(Th_cen+4);
        dy(Rot_cen+4) = F(Rot_cen)/(mass*la^2*100) - damp_u/(mass*la^2*10)*y(Rot_cen+4);
        if isnan(friction_tube)        
            fric_x(i,j) = frition_tube(i,j,1);
            fric_y(i,j) = frition_tube(i,j,2);
        end
        u_x(i,j) = y(U_cen);
        u_y(i,j) = y(V_cen);
        v_x(i,j) = y(U_cen+4);
        v_y(i,j) = y(V_cen+4);
        th(i,j) = y(Th_cen);
        rot(i,j) = y(Rot_cen);        
%         if i == 1 && j == 1 && seg == 1
%             progress
%         end  
    end
end
% y(Th_cen)
% friction_tube(1,:,1)
end

function dU = vector(theta,beta,edge)
global la
if strcmp(edge,'right')
    dUx = 2*la*cos(theta)*cos(beta);
    dUy = 2*la*cos(theta)*sin(beta);
elseif strcmp(edge,'left')
    dUx = -2*la*cos(theta)*cos(beta);
    dUy = -2*la*cos(theta)*sin(beta);
elseif strcmp(edge,'up')
    dUx = -2*la*cos(theta)*sin(beta);
    dUy = 2*la*cos(theta)*cos(beta);
elseif strcmp(edge,'down')
    dUx = 2*la*cos(theta)*sin(beta);
    dUy = -2*la*cos(theta)*cos(beta);
end
dU = [dUx;dUy];
end

function [F_input,M_input] = func_actuation(y,T,i,j)
    % basic information
    global Time_pre Time_switch Time_wait angle0
    global cycle actuator_tube Nx cycle_number
    k_M = 10;
    posi = ((i-1)*Nx + j-1)*8;
    
    F_input = 0; M_input = 0;
    [~,act_no] = max(i == actuator_tube(:,1) & j == actuator_tube(:,2));

    if T < Time_pre
        % no force applied, record the balacne angle
        M_input = 0;
        angle0(act_no) = y(posi+3);
    else
        % total time elapsed in this cycle
        cycle_time = mod(T-Time_pre,(Time_switch+Time_wait)*length(actuator_tube(:,1)));
        % cycle number in this cycle
        cycle_number = floor((T-Time_pre)/(Time_switch+Time_wait)/length(actuator_tube(:,1)));
        if cycle_number < cycle
            if cycle_time <= (Time_switch+Time_wait)*(act_no-1) % before actuation of this unit
                angle = (-1)^(cycle_number)*angle0(act_no);
                
            elseif (cycle_time <= (Time_switch+Time_wait)*(act_no-1)+Time_switch) && ...
                    (cycle_time > (Time_switch+Time_wait)*(act_no-1))
                % during the actuation
                angle = (angle0(act_no) -2*angle0(act_no)*(cycle_time-(Time_switch+Time_wait)*(act_no-1))...
                    /Time_switch)*(-1)^(cycle_number);
            else
                % after actuation
                angle = -(-1)^(cycle_number)*angle0(act_no);
            end
            M_input = k_M*(angle-y(posi+3));
        end
    end
end

function F_friction = func_friction(y_fric,T,i,j)
% basic information
global Time_pre mass Nx actuator_tube mass_act cycle_number thst
global wheel_angle y_0 mu_0 mu_1 pin_support mu_pin
F_friction = [0,0,0];
velo = [y_fric(5),y_fric(6)];
omega = y_fric(7);
Rot_angle =  y_fric(4);
if T > Time_pre% && norm(velo)>5e-5
    
    % record the initial value of each cycle
    %% calculate friciton
    % decompose velocity to along wheel and perpendicular to wheel
    abs_angle = wheel_angle(i,j)+Rot_angle;
    v_r = velo(1)*cos(abs_angle) + velo(2)*sin(abs_angle);
    v_t = -velo(1)*sin(abs_angle) + velo(2)*cos(abs_angle);
    %         F_r =  -0.5*mass*9.8*(abs(fric_coeff(y_fric(3)))*tanh((v_r)./5e-5));
    %         F_t =  -0.5*mass*9.8*(abs(fric_coeff(y_fric(3)+pi/2))*tanh((v_t)./5e-5));
    F_r = - mu_0*mass*9.8*((abs(sin(y_fric(3)))*tanh((v_r)./2e-5))...
        + 0*v_r*abs(sin(y_fric(3))));
    F_t = - mu_0*mass*9.8*(((abs(sin(y_fric(3)+pi/2)))*tanh((v_t)./2e-5))...
        + 0*v_t*abs(sin(y_fric(3)+pi/2)));
    
    F_friction(1) = F_r*cos(abs_angle) - F_t*sin(abs_angle);
    F_friction(2) = F_r*sin(abs_angle) + F_t*cos(abs_angle);
    
    %% rotational resistance
    posi = ((i-1)*Nx + j-1)*8;
    
    if y_0(posi+3)*(-1)^cycle_number>0
        fric_sign = 1;
    else
        fric_sign = -1;
    end
    index = max((thst-abs(y_fric(3))),0)./(thst/20);
    r = 0.02;
    F_friction(3) =r*9.8*mass*(tanh(index))*...
        ((-1)^cycle_number*fric_sign*0.95*tanh(norm(omega)./1e-3))-(tanh(index))*...
        mu_1*0.02*mass*9.8*omega;
    if sum(i == actuator_tube(:,1) & j == actuator_tube(:,2))
        F_friction(1:2) = F_friction(1:2).*(mass+mass_act)/mass; % with the actuator, the unit is heavier
        F_friction(3) = 0;
    end
    
    %%  for pin supports
    if pin_support == 1
        if norm(velo) == 0
            norm_fric = [0,0];
        else
            norm_fric = -velo./norm(velo);
        end
        F_friction(1:2) = mass*9.8*mu_pin*norm_fric*tanh(norm(velo)./1e-3);
        if j == 1
            F_friction(1:2) = 1*mass*9.8*mu_pin*norm_fric*tanh(norm(velo)./1e-3);
        end
        F_friction(3) = -0.005*omega;
    end
end
end