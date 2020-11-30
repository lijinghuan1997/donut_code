%% params_input
params.c=2.99792458e8;
params.e=1.60218e-19;
params.mp=1.6726231e-27;
params.me=9.1093897e-31;
params.m0=params.me/params.mp;
params.bz = 23;  
params.n0e = 87;
params.k=20;
params.delta = 0.4;
params.omega=0;
params.lambda=1.5;
params.te_=20;
params.tau_=205/params.te_;
params.tau=(235-params.delta*params.tau_*params.te_)/(1-params.delta)/params.te_/params.lambda;
params.te = params.te_*params.lambda;
params.ti_ = params.te_*params.tau_;
params.ti = params.te*params.tau;
params.z=1;
params.m = params.m0/params.z;
params.T0 = (params.te_*params.lambda*(params.z+params.tau))*params.e*1e7;
params.B0 = 1e5*sqrt(4*pi*params.n0e*params.T0); 
params.be=6.0d0/params.B0;
%params.be_=-5.8874/params.B0;
params.be_=0/params.B0;
params.n0i = params.n0e*(1-params.delta)/(1-params.be*params.B0/params.bz) +...
    params.n0e*params.delta/(1-params.be_*params.B0/params.bz);
params.di = sqrt(2.998e10^2*params.mp*1e3/(4*pi*params.n0e*(2.998e9*1.6e-19)^2));
params.D=params.di*params.k;
params.Omega0=-1/sqrt((params.D)^4/(3*10^10)^2/(params.T0)*4*pi*(params.e*3*10^9)^2*params.n0e);
params.Omegai = params.Omega0 * params.omega/(params.omega-1);
params.Omegae = params.Omega0/(params.omega-1);
params.ae=params.m/(params.m + params.omega^2);
params.ai=params.omega^2/(params.m + params.omega^2);
params.epsilon = (params.omega^2+params.m)*params.z/params.k^2/(params.omega-1)^2;
v0 = sqrt(params.te*params.e*2/params.me);
v0i = sqrt(params.ti*params.e*2/params.mp);
%% load 初始电磁场
%inputs = [0:0.001:5 5.5:0.5:500 550:50:10000];
%inputs=[0:0.01:0.5];
%inputs=[0:0.001:0.2,0.21:0.01:5 5.5:0.5:500 550:50:10000];
%inputs=[0:0.01:5 5.5:0.5:500 550:50:10000];
inputs=[0:0.01:0.5,0.51:0.01:6,6.1:0.1:10];
sol = ode45(@(x,y)vpt2(x,y,params),inputs, [0;-params.bz/params.B0;0]);
[y, yp] = deval(inputs,sol);
x = sqrt(2*inputs);
E = -yp(3,:).*x*(params.te_*params.lambda+params.te_*params.lambda*params.tau)/params.D*100;
B = -y(2,:)*params.B0;
data={y,x,E,B};
disp('initial state finished')
%data={zeros(3,length(inputs)),x,zeros(1,length(inputs)),17*ones(1,length(inputs))};
%% ode+电磁场存储
jieshumax=2;
for i=1:jieshumax
    %inputs = [0:0.00001:5 5.5:0.5:500 550:50:10000];
    %inputs=[0:0.01:0.5];
    a=max(data{i,2});
    a=a-0.9;
    a=a^2/2;
%plot    
    %inputs=[0:0.001:0.2,0.21:0.01:5 5.5:0.5:500 550:50:a];
    %inputs=[0:0.01:5 5.5:0.5:500 550:50:a];
    %inputs=[0:0.01:2];
%测试    
    inputs=[0:0.001:0.5,0.51:0.01:6,6.1:0.1:a];
    sol = ode45(@(xx,y)vptl3(xx,y,data{i,2},data{i,3},data{i,4},params),inputs, [0;-params.bz/params.B0;0]);
    [y, yp] = deval(inputs,sol);
    x = sqrt(2*inputs);
    E = -yp(3,:).*x*(params.te_*params.lambda+params.te_*params.lambda*params.tau)/params.D*100;
    B = -y(2,:)*params.B0;
    data{i+1,1}=y;
    data{i+1,2}=x;
    data{i+1,3}=E;
    data{i+1,4}=B;
    disp(i)
end
disp('ode finished')
%plot(data{2,2},data{2,4});
jieshu=jieshumax;
x=data{jieshu,2};
E=data{jieshu,3};
B=data{jieshu,4};
y=data{jieshu+1,1};
xnew=data{jieshu+1,2};
Enew=data{jieshu+1,3};
Bnew=data{jieshu+1,4};
rr=@(xx,yy)sqrt(xx.^2+yy.^2);
E_e=@(xx)interp1(x,E,xx,'spline');
B_e=@(xx)interp1(x,B,xx,'spline');
        rc=@(vx,vy,xx,yy)params.me*v0/B_e(rr(xx,yy))*1e9/params.e/params.D*100*sqrt(vx.^2+vy.^2); %guiyihua
        rc_x=@(vx,vy,xx,yy)xx+rc(vx,vy,xx,yy).*(-vy).*judge_0(vx,vy);
        rc_y=@(vx,vy,xx,yy)yy+rc(vx,vy,xx,yy).*vx.*judge_0(vx,vy);
        r_=@(vx,vy,xx,yy)sqrt(rc_x(vx,vy,xx,yy).^2+rc_y(vx,vy,xx,yy).^2);
        % function handles to calcute miu
        vE=@(vx,vy,xx,yy)-E_e(r_(vx,vy,xx,yy))./B_e(r_(vx,vy,xx,yy))*1e9;
        vf=@(vx,vy,vz,xx,yy)(vx*(-yy)+vy*xx).*judge_0(xx,yy);
        B_tidu=dif2(B,x)/params.D*100*1e-9;
        B_tidu_e=@(r)interp1(x,B_tidu,r,'spline');
        vD=@(vx,vy,xx,yy)1/2*params.me*(vx.^2+vy.^2)*v0^2./params.e./B_e(r_(vx,vy,xx,yy)).^2 ...
            .*B_tidu_e(r_(vx,vy,xx,yy))*1e18;
        vdrift=@(vx,vy,xx,yy)vE(vx,vy,xx,yy)-vD(vx,vy,xx,yy);
        vdrift_x=@(vx,vy,xx,yy)vdrift(vx,vy,xx,yy).*(-rc_y(vx,vy,xx,yy)).*judge_0(rc_x(vx,vy,xx,yy),rc_y(vx,vy,xx,yy));
        vdrift_y=@(vx,vy,xx,yy)vdrift(vx,vy,xx,yy).*(rc_x(vx,vy,xx,yy)).*judge_0(rc_x(vx,vy,xx,yy),rc_y(vx,vy,xx,yy));
miu_guiyi1=@(vx,vy,xx,yy)((vx*v0-vdrift_x(vx,vy,xx,yy)).^2+(vy*v0-vdrift_y(vx,vy,xx,yy)).^2) ...
    *params.be*params.B0./interp1(xnew,Bnew,rr(xx,yy),'spline')/v0^2;
miu_guiyi2=@(vx,vy,xx,yy)((vx*v0-vdrift_x(vx,vy,xx,yy)).^2+(vy*v0-vdrift_y(vx,vy,xx,yy)).^2) ...
    *params.be_*params.B0./interp1(xnew,Bnew,rr(xx,yy),'spline')/v0^2;
vf=@(vx,vy,vz,xx,yy)(vx*(-yy)+vy*xx).*judge_0(xx,yy);
f_e_rec_1= @(vx,vy,xx,yy)params.n0e*(1-params.delta)*(1/pi)*exp(...
        -vx.^2 - vy.^2 + miu_guiyi1(vx,vy,xx,yy) ...
        + sign(params.Omegae)*vf(vx,vy,0,xx,yy).*rr(xx,yy)*sqrt(2*params.ae*params.epsilon*(1+params.tau)) +...
        (1+params.tau)*(interp1(xnew,y(3,:),rr(xx,yy),'spline')-1/(params.omega-1)*interp1(xnew,y(1,:),rr(xx,yy),'spline'))); 
f_e_rec_2=@(vx,vy,xx,yy)params.n0e*params.delta*(params.lambda/pi)*exp(...
        -params.lambda*(vx.^2+vy.^2-miu_guiyi2(vx,vy,xx,yy)) ...
        + params.lambda*(1+params.tau)*interp1(xnew,y(3,:),rr(xx,yy),'spline'));
f_e_rec__=@(vx,vy,xx,yy)f_e_rec_1(vx,vy,xx,yy)+f_e_rec_2(vx,vy,xx,yy);
nn_e=@(xx,yy)integral2(@(vx,vy)f_e_rec__(vx,vy,xx,yy),-3,3,-3,3);
params.n0i=nn_e(0,0);
f_i_rec__=@(vx,vy,xx,yy)params.n0i*(1-params.delta)*(1/pi)^1*exp(...
        -(vx.^2+vy.^2)+sign(params.Omegai)*vf(vx,vy,0,xx,yy).*rr(xx,yy)*sqrt(2*params.ai*params.epsilon*(1+params.tau)/params.tau)+...
        (1+params.tau)/params.tau*(params.omega/(params.omega-1)*interp1(xnew,y(1,:),rr(xx,yy),'spline')-interp1(xnew,y(3,:),rr(xx,yy),'spline')))+...
            params.n0i*params.delta*(params.lambda*params.tau/pi/params.tau_)^1*exp(...
        -params.lambda*params.tau/params.tau_*(vx.^2+vy.^2)-...
        params.lambda/params.tau_*(1+params.tau)*interp1(xnew,y(3,:),rr(xx,yy),'spline')); 
    
nn_i=@(xx,yy)integral2(@(vx,vy)f_i_rec__(vx,vy,xx,yy),-3,3,-3,3); 
g_e = @(v,t,p,xx)1e-4*1e6*v0^1*v.^4/2.*f_e(v.*sin(t).*cos(p),...
    v.*sin(t).*sin(p), v.*cos(t), xx);
averaged_g = @(v,t,xx)integral(@(p)g_e(v,t,p,xx),0,2*pi)/(2*pi);
g_i = @(v,t,p,xx)1e-4*1e6*v0i^1*v.^4/2.*f_i(v.*sin(t).*cos(p),...
    v.*sin(t).*sin(p), v.*cos(t), xx);
averaged_gi = @(v,t,xx)integral(@(p)g_i(v,t,p,xx),0,2*pi)/(2*pi);
%% save data to fortran
Breal=interp1(xnew,Bnew,[0:1:500000]/params.D*10^2,'spline');
Ereal=interp1(xnew,Enew,[0:1:500000]/params.D*10^2,'spline');
Ychazhi3=interp1(xnew,y(3,:),[0:1:500000]/params.D*10^2,'spline');
Ychazhi1=interp1(xnew,y(1,:),[0:1:500000]/params.D*10^2,'spline');
Ychazhi2=interp1(xnew,y(2,:),[0:1:500000]/params.D*10^2,'spline');
Rreal=0:1:500000;
A=[Rreal',Breal',Ereal',Ychazhi1',Ychazhi2',Ychazhi3'];
save data_miu_diedai.txt -ascii A
%% function    
function dydx=vptl3(xx,y,x,E,B,params)
        dydx=zeros(3,1);  
        lambda = params.lambda;
        tau = params.tau;
        tau_ = params.tau_;
        delta = params.delta;
        k = params.k;
        m = params.m/1;
        z = params.z;
        be = params.be;
        be_ = params.be_;
        qujian=3;
        v0 = sqrt(params.te*1.6e-19*2/params.me);
        rr=@(xx,yy)sqrt(xx.^2+yy.^2);
        E_e=@(xx)interp1(x,E,xx,'spline');
        B_e=@(xx)interp1(x,B,xx,'spline');
        rc=@(vx,vy,xx,yy)params.me*v0/B_e(rr(xx,yy))*1e9/params.e/params.D*100*sqrt(vx.^2+vy.^2); %guiyihua
        rc_x=@(vx,vy,xx,yy)xx+rc(vx,vy,xx,yy).*(-vy).*judge_0(vx,vy);
        rc_y=@(vx,vy,xx,yy)yy+rc(vx,vy,xx,yy).*vx.*judge_0(vx,vy);
        r_=@(vx,vy,xx,yy)sqrt(rc_x(vx,vy,xx,yy).^2+rc_y(vx,vy,xx,yy).^2);
        % function handles to calcute miu
        vE=@(vx,vy,xx,yy)-E_e(r_(vx,vy,xx,yy))./B_e(r_(vx,vy,xx,yy))*1e9;
        vf=@(vx,vy,vz,xx,yy)(vx*(-yy)+vy*xx).*judge_0(xx,yy);
        B_tidu=dif2(B,x)/params.D*100*1e-9;
        B_tidu_e=@(r)interp1(x,B_tidu,r,'spline');
        vD=@(vx,vy,xx,yy)1/2*params.me.*(vx.^2+vy.^2)*v0^2./params.e./B_e(r_(vx,vy,xx,yy)).^2 ...
            .*B_tidu_e(r_(vx,vy,xx,yy))*1e18; %./sqrt(1-(vx.^2+vy.^2).*v0^2./params.c.^2)
        vdrift=@(vx,vy,xx,yy)vE(vx,vy,xx,yy)-vD(vx,vy,xx,yy);
        vdrift_x=@(vx,vy,xx,yy)vdrift(vx,vy,xx,yy).*(-rc_y(vx,vy,xx,yy)).*judge_0(rc_x(vx,vy,xx,yy),rc_y(vx,vy,xx,yy));
        vdrift_y=@(vx,vy,xx,yy)vdrift(vx,vy,xx,yy).*(rc_x(vx,vy,xx,yy)).*judge_0(rc_x(vx,vy,xx,yy),rc_y(vx,vy,xx,yy));
        miu_guiyi1=@(vx,vy,xx,yy,y2)((vx*v0-vdrift_x(vx,vy,xx,yy)).^2+(vy*v0-vdrift_y(vx,vy,xx,yy)).^2) ...
                *params.be./(-y2)/v0^2;
        miu_guiyi2=@(vx,vy,xx,yy,y2)((vx*v0-vdrift_x(vx,vy,xx,yy)).^2+(vy*v0-vdrift_y(vx,vy,xx,yy)).^2) ...
                *params.be_./(-y2)/v0^2;
        f_e_rec_1= @(vx,vy,xx,yy,y2,y1,y3)params.n0e*(1-params.delta)*(1/pi)*exp(...
            -vx.^2 - vy.^2 + miu_guiyi1(vx,vy,xx,yy,y2) ...
            + sign(params.Omegae)*vf(vx,vy,0,xx,yy).*rr(xx,yy)*sqrt(2*params.ae*params.epsilon*(1+params.tau)) +...
            (1+params.tau)*(y3-1/(params.omega-1)*y1)); 
        f_e_rec_2= @(vx,vy,xx,yy,y2,y1,y3)params.n0e*params.delta*(params.lambda/pi)*exp(...
            -params.lambda*(vx.^2+vy.^2-miu_guiyi2(vx,vy,xx,yy,y2)) ...
            + params.lambda*(1+params.tau)*y3);
        f_e_rec__=@(vx,vy,xx,yy,y2,y1,y3)f_e_rec_1(vx,vy,xx,yy,y2,y1,y3)+f_e_rec_2(vx,vy,xx,yy,y2,y1,y3);
    Ee=integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*xx),0,y(2),y(1),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
    params.n0i=integral2(@(vx,vy)f_e_rec__(vx,vy,1e-7,0,-params.bz/params.B0,0,0),-qujian,qujian,-qujian,qujian);
    nie = params.n0i/params.n0e;
    if Ee==inf && xx>0.5
        while Ee>2
        qujian=qujian-0.1;
        Ee=integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*xx),0,y(2),y(1),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
        end
    end
    Ee_func=@(t)integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*t),0,y(2),y(1),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
    Ee_func_y2=@(t)integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*xx),0,t,y(1),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
    Ee_func_y1=@(t)integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*xx),0,y(2),t,y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
    Ee_func_y3=@(t)integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*xx),0,y(2),y(1),t),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
    Ee_=integral2(@(vx,vy)f_e_rec_2(vx,vy,sqrt(2*xx),0,y(2),y(1),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(params.delta);
    Ee__func=@(t)integral2(@(vx,vy)f_e_rec_2(vx,vy,sqrt(2*t),0,y(2),y(1),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(params.delta);
    Ee__func_y3=@(t)integral2(@(vx,vy)f_e_rec_2(vx,vy,sqrt(2*xx),0,y(2),y(1),t),-qujian,qujian,-qujian,qujian)/params.n0e/(params.delta);
    Ee__func_y2=@(t)integral2(@(vx,vy)f_e_rec_2(vx,vy,sqrt(2*xx),0,t,y(1),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(params.delta);
    Ei=exp((z+tau)/tau*(params.ai*params.epsilon*xx+(params.omega/(params.omega-1))*y(1)-y(3)));
    Ei_=exp(-lambda*(1+tau)/tau_*y(3));
    dydx(1)=y(2);
    if xx==0
        xx_=1e-7;
        je =integral2(@(vx,vy)f_e_rec__(vx,vy,sqrt(2*xx_),0,y(2),y(1),y(3)).*vy*v0,-qujian,qujian,-qujian,qujian);
        
        dydx(2)=(1-delta)*(je/params.Omega0/sqrt(xx_*params.D^2*2/1e4)/params.n0e/(1-delta) - nie*Ei*params.omega/(params.omega-1));
        Psi=(1+tau)/(params.omega-1)*(nie*params.omega/tau*Ei)-(Ee_func_y1(y(1)+1e-7)-Ee)/1e-7;
        PPsi=-(Ee_func_y2(y(2)+1e-7)-Ee)/1e-7-(Ee__func_y2(y(2)+1e-7)-Ee_)/1e-7*delta/(1-delta);
        Theta=params.epsilon*(1+tau)*(nie*params.ai*Ei/tau)-(Ee_func(xx+1e-7)-Ee)/1e-7-(Ee__func(xx+1e-7)-Ee_)/1e-7*delta/(1-delta);   
        K=(Ee_func_y3(y(3)+1e-7)-Ee)/1e-7+(Ee__func_y3(y(3)+1e-7)-Ee_)/1e-7*delta/(1-delta)+nie/tau*Ei*(1+tau)+(1+tau)*nie*lambda/tau_*delta/(1-delta)*Ei_;  
        dydx(3) = (Psi*y(2)+PPsi*dydx(2)+Theta)/K;
    else
        
        je =integral2(@(vx,vy)f_e_rec__(vx,vy,sqrt(2*xx),0,y(2),y(1),y(3)).*vy*v0,-qujian,qujian,-qujian,qujian);
        dydx(2)=(1-delta)*(je/params.Omega0/sqrt(xx*params.D^2*2/1e4)/params.n0e/(1-delta) - nie*Ei*params.omega/(params.omega-1));
        while abs(dydx(2))>1
            je =integral2(@(vx,vy)f_e_rec__(vx,vy,sqrt(2*xx),0,y(2),y(1),y(3)).*vy*v0,-qujian,qujian,-qujian,qujian);
            dydx(2)=(1-delta)*(je/params.Omega0/sqrt(xx*params.D^2*2/1e4)/params.n0e/(1-delta) - nie*Ei*params.omega/(params.omega-1));
            qujian=qujian-0.1;
        end
        Psi=(1+tau)/(params.omega-1)*(nie*params.omega/tau*Ei)-(Ee_func_y1(y(1)+1e-7)-Ee)/1e-7;
        PPsi=-(Ee_func_y2(y(2)+1e-7)-Ee)/1e-7-(Ee__func_y2(y(2)+1e-7)-Ee_)/1e-7*delta/(1-delta);
        Theta=params.epsilon*(1+tau)*(nie*params.ai*Ei/tau)-(Ee_func(xx+1e-7)-Ee)/1e-7-(Ee__func(xx+1e-7)-Ee_)/1e-7*delta/(1-delta);   
        K=(Ee_func_y3(y(3)+1e-7)-Ee)/1e-7+(Ee__func_y3(y(3)+1e-7)-Ee_)/1e-7*delta/(1-delta)+nie/tau*Ei*(1+tau)+(1+tau)*nie*lambda/tau_*delta/(1-delta)*Ei_;  
        dydx(3) = (Psi*y(2)+PPsi*dydx(2)+Theta)/K;
    end    
end
function dydx=vpt2(x,y,params)
    % params: 常数结构体，含有下列参数
    lambda = params.lambda;
    tau = params.tau;
    tau_ = params.tau_;
    delta = params.delta;
    k = params.k;
    m = params.m/1;
    z = params.z;
    be = params.be;
    be_ = params.be_;
    params.n0i = params.n0e*(1-params.delta)/(1-params.be*params.B0/params.bz) +...
    params.n0e*params.delta/(1-params.be_*params.B0/params.bz);
    nie = params.n0i/params.n0e;
    
    dydx=zeros(3,1);
    
%     k1=delta*lambda/(1-delta);
%     k2=z/tau;
%     k3=z*lambda/tau_*delta/(1-delta);
%     m1=lambda*(z+tau);
%     m2=-z*(z+tau)/tau;
%     m3=-lambda*(z+tau)/tau_;
    xie=be/y(2);
    xie_ = be_/y(2);
    Ee=exp((z+tau)*(params.ae*params.epsilon*x/(1+xie)-(1/(params.omega-1))*y(1)+y(3)));
    Ei=exp((z+tau)/tau*(params.ai*params.epsilon*x+(params.omega/(params.omega-1))*y(1)-y(3)));
    Ee_=exp(lambda*(1+tau)*y(3));
    Ei_=exp(-lambda*(1+tau)/tau_*y(3));
    
    dydx(1)=y(2);
    dydx(2)=(1-delta)*(Ee/(params.omega-1)/(1+xie)^2 - nie*Ei*params.omega/(params.omega-1));
    % assume: be is small enough,
    % dydx(3)=-1/(z+tau)*...
        % (-(z+tau)*Ee*y(2)+z*m/k^2*(z+tau)*Ee/(1+xie)+...
        % (1-delta)*(exp(m1*y(3))-Ee-Ee*z*m/k^2*x*(z+tau)/(1+xie))*dydx(2))/... 
        % (Ee+k1*exp(m1*y(3))+k2*exp(m2*y(3))+k3*exp(m3*y(3))); % K;
    % assume: bi=0
    
    % TODO: equation 3
    Psi=(1+tau)/(params.omega-1)*(nie*params.omega/tau*Ei+Ee/(1+xie));
    PPsi=-Ee*(1+(1+tau)*params.ae*params.epsilon*x/(1+xie))*(be/(be+y(2))^2)...
        - Ee_*delta/(1-delta)*(be_/(be_+y(2))^2);
    Theta=params.epsilon*(1+tau)*(nie*params.ai*Ei/tau-params.ae*Ee/(1+xie)^2);
    K=(1+tau)*(Ee/(1+xie)+lambda*delta/(1-delta)/(1+xie_)*Ee_+nie/tau*Ei+nie*lambda/tau_*delta/(1-delta)*Ei_);
    
%     dydx(3) = (Ee*y(2)/(1+xie) - m/k^2*Ee/(1+xie)^2 - Ee*(1/(1+tau)+m/k^2*x/(1+xie))*(be/(be+y(2))^2)*dydx(2)) /...
%         (nie*(k2*exp(m2*y(3)) + k3*exp(m3*y(3))) + (k1*exp(m1*y(3)) + Ee/(1+xie)));
    dydx(3) = (Psi*y(2)+PPsi*dydx(2)+Theta)/K;
end

function result=judge_0(x,y)
    if x==0 & y==0
        result=0;
    else
        result=1./sqrt(x.^2+y.^2);
    end
end