% ***************************************************************************
%  Time advancement code
%   Steps a single timestep for u, v, c, rho, q2, q2l, l, kz, nu_t, kq
%   All diffusion/viscous terms handled implicitly
% ***************************************************************************

% ***************************************************************************
%  Update pressure forcing term for the current timestep
% ***************************************************************************
for i = 1:N
   if T_Px == 0.0
      Px(i) = Px0;  % Steady and constant forcing for now
   else
      Px(i) = Px0*cos(2*pi*t(m)/(3600.*T_Px));
   end
end

%Update shear velocity at bottom boundary for use later
ustar = abs(U(1))*sqrt(C_D); % Explicit dependence on C_D

%****************************************************************************
% Update parameters for the model, Sm and Sh
%****************************************************************************

for i=1:N		
   Gh=-(N_BV(i)*L(i)/(Q(i)+SMALL))^2; 
   %set LIMITER for Gh 
   Gh=min(Gh, 0.0233);
   Gh=max(Gh, -0.28);
   num=B1^(-1/3)-A1*A2*Gh*((B2-3*A2)*(1-6*A1/B1)-3*C1*(B2+6*A1));
   dem=(1-3*A2*Gh*(B2+6*A1))*(1-9*A1*A2*Gh);
   Sm(i)=num/dem;
   Sh(i)=A2*(1-6*A1/B1)/(1-3*A2*Gh*(B2+6*A1)); 
end

% ***************************************************************************
%  Place previous variable f into fp (i.e. u into up, q2 into q2p, etc)
% ***************************************************************************
Up = U;
Vp = V;
Cp = C;
Q2p = Q2;
Q2Lp = Q2L;
Lp = L;
Kzp = Kz;
nu_tp = nu_t;
Kqp = Kq;
N_BVp = N_BV;

% ***************************************************************************
%  Advance velocity (U,V)
% ***************************************************************************
for i = 2:N-1
   aU(i) = -0.5*beta*(nu_tp(i)+nu_tp(i-1));
   bU(i) = 1+0.5*beta*(nu_tp(i+1)+2*nu_tp(i)+nu_tp(i-1));
   cU(i) = -0.5*beta*(nu_tp(i)+nu_tp(i+1));
   dU(i) = Up(i) - dt*Px(i);
end
%bottom boundary: log-law
bU(1) = 1+0.5*beta*(nu_tp(2)+nu_tp(1)+2*(sqrt(C_D)/kappa)*nu_tp(1));
cU(1) = -0.5*beta*(nu_tp(2)+nu_tp(1));
dU(1) = Up(1) - dt*Px(1);
%top boundary: no stress
aU(N) = -0.5*beta*(nu_tp(N)+nu_tp(N-1));
bU(N) = 1+0.5*beta*(nu_tp(N)+nu_tp(N-1));
dU(N) = Up(N) - dt*Px(N);
% Use Thomas algorithm to solve for U
for i=2:N
   bU(i) = bU(i) - aU(i)/bU(i-1)*cU(i-1);
   dU(i) = dU(i) - aU(i)/bU(i-1)*dU(i-1);
end
U(N) = dU(N)/bU(N);
for i=N-1:-1:1
  U(i) = 1/bU(i)*(dU(i)-cU(i)*U(i+1));
end


% ***************************************************************************
%  Advance scalars/density (C, rho) 
% ***************************************************************************
for i = 2:N-1
   aC(i) = -0.5*beta*(Kzp(i)+Kzp(i-1));
   bC(i) = 1+0.5*beta*(Kzp(i+1)+2*Kzp(i)+Kzp(i-1));
   cC(i) = -0.5*beta*(Kzp(i)+Kzp(i+1));
   dC(i) = Cp(i);
end
%bottom boundary: no flux for scalars
bC(1) = 1+0.5*beta*(Kzp(2)+Kzp(1));
cC(1) = -0.5*beta*(Kzp(2)+Kzp(1));
dC(1) = Cp(1);
%top boundary: no flux for scalars
aC(N) = -0.5*beta*(Kzp(N)+Kzp(N-1));
bC(N) = 1+0.5*beta*(Kzp(N)+Kzp(N-1));
dC(N) = Cp(N);
% Use Thomas algorithm to solve for C
for i=2:N
   bC(i) = bC(i) - aC(i)/bC(i-1)*cC(i-1);
   dC(i) = dC(i) - aC(i)/bC(i-1)*dC(i-1);
end
C(N) = dC(N)/bC(N);
for i=N-1:-1:1
  C(i) = 1/bC(i)*(dC(i)-cC(i)*C(i+1));
end
%update density and Brunt-Vaisala frequency
for i = 1:N
   rho(i)=rho0*(1-alpha*(C(i)-15)); % Single scalar, linear equation of state
end
N_BV(1)=sqrt((-g/rho0)*(rho(2)-rho(1))/(dz));
for i =2:N-1
   N_BV(i)=sqrt((-g/rho0)*(rho(i+1)-rho(i-1))/(2*dz));
end
N_BV(N)=sqrt((-g/rho0)*(rho(N)-rho(N-1))/(dz));



% ***************************************************************************
%  Advance turbulence parameters (q2, q2l - q2 first, then q2l)
% ***************************************************************************
diss = np.zeros(N)
for i=2:N-1   
   diss = 2*dt*(Q2p(i))^(1/2)/(B1*Lp(i)); % coefficient for linearized term
   aQ2(i) = -0.5*beta*(Kqp(i)+Kqp(i-1));
   bQ2(i) = 1+0.5*beta*(Kqp(i+1)+2*Kqp(i)+Kqp(i-1))+diss;
   cQ2(i) = -0.5*beta*(Kqp(i)+Kqp(i+1));
   dQ2(i) = Q2p(i) + 0.25*beta*(nu_tp(i))*(Up(i+1)-Up(i-1))^2 - dt*(Kzp(i))*(N_BVp(i))^2;
end
%Bottom boundary (i=1)
Q2bot=B1^(2/3)*ustar^2; %Q2(0) is prescribed
bndryterm = 0.5*beta*(Kqp(1))*Q2bot;
diss = 2*dt*(Q2p(1))^(1/2)/(B1*Lp(1));
bQ2(1) = 1+0.5*beta*(Kqp(2)+Kqp(1))+diss;
cQ2(1) = -0.5*beta*(Kqp(1)+Kqp(2));
dQ2(1) = Q2p(1) + dt*(ustar^4)/nu_tp(i)- dt*(Kzp(1))*(N_BVp(1))^2 + bndryterm;
%Top boundary (i=N)
diss = 2*dt*(Q2p(N))^(1/2)/(B1*Lp(N));
aQ2(N) = -0.5*beta*(Kqp(N)+Kqp(N-1));
bQ2(N) = 1+0.5*beta*(Kqp(N)+2*Kqp(N)+Kq(N-1))+diss;
dQ2(N) = Q2p(N) + 0.25*beta*(nu_tp(N))*(Up(N)-Up(N-1))^2 - 4*dt*(Kzp(N))*(N_BVp(N))^2 ;
% Use Thomas algorithm to solve for Q2
for i=2:N
   bQ2(i) = bQ2(i) - aQ2(i)/bQ2(i-1)*cQ2(i-1);
   dQ2(i) = dQ2(i) - aQ2(i)/bQ2(i-1)*dQ2(i-1);
end
Q2(N) = dQ2(N)/bQ2(N);
for i=N-1:-1:1
	Q2(i) = 1/bQ2(i)*(dQ2(i)-cQ2(i)*Q2(i+1));
end

%  Kluge to prevent negative values from causing instabilities
for i=1:N
   if Q2(i)<SMALL
      Q2(i)=SMALL;
   end
end
% *******************************************************************

for i=2:N-1
   diss=2*dt*(Q2p(i))^(1/2)/(B1*Lp(i))*(1+E2*(Lp(i)/(kappa*abs(-H-z(i))))^2+E3*(Lp(i)/(kappa*abs(z(i))))^2);
   aQ2L(i) = -0.5*beta*(Kqp(i)+Kqp(i-1));
   bQ2L(i) = 1+0.5*beta*(Kqp(i+1)+2*Kqp(i)+Kqp(i-1))+diss;
   cQ2L(i) = -0.5*beta*(Kqp(i)+Kqp(i+1));
   dQ2L(i) = Q2Lp(i) + 0.25*beta*(nu_tp(i))*E1*Lp(i)*(Up(i+1)-Up(i-1))^2 ...
		- 2*dt*Lp(i)*E1*(Kzp(i))*(N_BVp(i))^2;
end
%Bottom boundary (i=1)
Q2Lbot= B1^(2/3)*ustar^2*kappa*zb; 
bndryterm= 0.5*beta*(Kqp(1))*Q2Lbot;
diss=2*dt*(Q2p(1))^(1/2)/(B1*Lp(1))*(1+E2*(Lp(1)/(kappa*abs(-H-z(1))))^2+E3*(Lp(1)/(kappa*abs(z(1))))^2);
bQ2L(1) = 1+0.5*beta*(Kqp(2)+Kqp(1))+diss;
cQ2L(1) = -0.5*beta*(Kq(1)+Kq(2));
dQ2L(1) = Q2Lp(1) + dt*((ustar^4)/nu_tp(1))*E1*Lp(1)-dt*Lp(1)*E1*(Kzp(1))*(N_BVp(1))^2+bndryterm;
%Top boundary (i=N)
diss=2*dt*(Q2p(N))^(1/2)/(B1*Lp(N))*(1+E2*(Lp(N)/(kappa*abs(-H-z(N))))^2+E3*(Lp(N)/(kappa*abs(z(N))))^2);
aQ2L(N) = -0.5*beta*(Kqp(N)+Kqp(N-1));
bQ2L(N) = 1+0.5*beta*(Kqp(N)+2*Kqp(N)+Kqp(N-1))+diss;
dQ2L(N) = Q2Lp(N) + 0.25*beta*(nu_tp(N))*E1*Lp(N)*(Up(N)-Up(N-1))^2 - 2*dt*Lp(N)*E1*(Kzp(N))*(N_BVp(N))^2;
% Use Thomas algorithm to solve for Q2L
for i=2:N
   bQ2L(i) = bQ2L(i) - aQ2L(i)/bQ2L(i-1)*cQ2L(i-1);
   dQ2L(i) = dQ2L(i) - aQ2L(i)/bQ2L(i-1)*dQ2L(i-1);
end
Q2L(N) = dQ2L(N)/bQ2L(N);
for i=N-1:-1:1
   Q2L(i) = 1/bQ2L(i)*(dQ2L(i)-cQ2L(i)*Q2L(i+1));
end

%  Kluge to prevent negative values from causing instabilities
for i=1:N
   if Q2L(i)<SMALL
      Q2L(i)=SMALL;
   end	
end

% ***************************************************************************
%  Calculate turbulent lengthscale (l) and mixing coefficients (kz, nu_t, kq)
%     Works will all updated values 
% ***************************************************************************

for i=1:N		
   Q(i)=sqrt(Q2(i));
   L(i)=Q2L(i)/(Q2(i)+SMALL); 
   %limit due to stable stratification
   if (L(i)^2*(N_BV(i))^2>0.281*Q2(i)) 
      	%Adjust Q2L as well as L
      	Q2L(i)=Q2(i)*sqrt(0.281*Q2(i)/((N_BV(i))^2+SMALL));
      	L(i)=Q2L(i)/Q2(i);
   end
   %Keep L from becoming zero -- zb=bottom roughness parameter
   if  (abs(L(i))<=zb)
       	L(i)=zb;
   end
   %update diffusivities
   Kq(i)=Sq*Q(i)*L(i) + nu; 
   nu_t(i)=Sm(i)*Q(i)*L(i) + nu; 
   Kz(i)=Sh(i)*Q(i)*L(i) + nu;
end


