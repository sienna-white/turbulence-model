from numba import jit 
import numpy as np
import math
from constants import * 


@jit 
def TDMA(aX, bX, cX, dX, N):
    # Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
    # a = Lower Diag, b = Main Diag, c = Upper Diag, d = solution vector
    x = np.zeros(N)
    for i in range(1, N):
        bX[i] = bX[i] - aX[i]/bX[i-1]*cX[i-1]
        dX[i] = dX[i] - aX[i]/bX[i-1]*dX[i-1]
    x[-1] = dX[-1]/bX[-1]
    for i in range(N-2, -1, -1):
        x[i] = (1/bX[i])*(dX[i] - cX[i]*x[i+1])
    return x

def initialize_abcd(N):
    a = np.zeros(N)
    b = np.zeros(N)
    c = np.zeros(N)
    d = np.zeros(N)
    return a, b, c, d

def advance_velocity(self, Up, nu_tp, Px, W=None):

    aU, bU, cU, dU = initialize_abcd(self.N) 
    top = self.top
    beta = self.beta 
    dt  = self.dt 

    aU[1:top] = -beta/2*(nu_tp[1:top] + nu_tp[0:top-1])
    bU[1:top] = 1 + beta/2*(nu_tp[2:top+1] + 2*nu_tp[1:top] + nu_tp[0:top-1])
    cU[1:top] = -beta/2*(nu_tp[1:top] + nu_tp[2:top+1])
    dU[1:top] = Up[1:top] - dt*Px[1:top]

    # Bottom boundary: log-law
    bU[0] = 1 + beta/2*(nu_tp[1] + nu_tp[0] + 2*(math.sqrt(C_D)/kappa)*nu_tp[0])
    cU[0] = -beta/2*(nu_tp[1] + nu_tp[0])
    dU[0] = Up[0] - dt*Px[0]

    aU[top] = -beta/2*(nu_tp[top]+nu_tp[top-1])
    bU[top] = 1 + beta/2*(nu_tp[top]+nu_tp[top-1])
    dU[top] = Up[top] - dt*Px[top] + W # beta*(nu_tp[top]/2)*W

    # # Top boundary: no stress
    # aU[top] = -beta/2*(nu_tp[top]+nu_tp[top-1])
    # bU[top] = 1 + beta/2*(nu_tp[top]+nu_tp[top-1])
    # dU[top] = Up[top] - dt*Px[top]
    U = TDMA(aU, bU, cU, dU, self.N)
    return U


def advance_algae(self, ws, gamma, Kzp, Ap):
    aA, bA, cA, dA = initialize_abcd(self.N)
    wsdtdz = abs(ws* self.dt)/self.dz
    top = self.top
    beta = self.beta 
    dt  = self.dt 

    # If settling speed is UPWARD (swimming!)
    if ws>0:
        aA[1:top]  = -wsdtdz - beta/2 * (Kzp[0:top-1]+ Kzp[1:top]) 
        bA[1:top]  = 1 + wsdtdz - gamma[1:top]*dt  + beta/2*(Kzp[2:top+1] + 2*Kzp[1:top] + Kzp[0:top-1]) 
        cA[1:top]  = -beta/2 * (Kzp[1:top] + Kzp[2:top+1])
        dA = Ap

        # Bottom-Boundary: no flux for scalars
        bA[0] =  1 - (gamma[0]*dt) + beta/2*(Kzp[1] + Kzp[0]) + wsdtdz
        cA[0] = -beta/2 * (Kzp[1] + Kzp[0])
        dA[0] =  Ap[0]

        # Top-Boundary: no flux for scalars
        aA[top] = -wsdtdz - beta/2 * (Kzp[top] + Kzp[top-1])
        bA[top] = 1  - gamma[top]*dt + beta/2 * (Kzp[top] + Kzp[top-1]) # removed + wsdtdz 
        dA[top] = Ap[top]

    # If settling speed is DOWNWARD (sinking!)
    if ws<=0:
        aA[1:top]  = -beta/2 * (Kzp[0:top-1] + Kzp[1:top]) 
        bA[1:top]  = 1 + wsdtdz - gamma[1:top]*dt  + beta/2*(Kzp[2:top+1] + 2*Kzp[1:top] + Kzp[0:top-1]) 
        cA[1:top]  = -wsdtdz - beta/2 * (Kzp[1:top] + Kzp[2:top+1])
        dA = Ap

        # Bottom-Boundary: no flux for scalars
        bA[0] =  1 + wsdtdz - (gamma[0]*dt) + beta/2*(Kzp[1] + Kzp[0]) 
        cA[0] =  -wsdtdz -beta/2 * (Kzp[1] + Kzp[0])
        dA[0] =  Ap[0]

        # Top-Boundary: no flux for scalars
        aA[top] =  -beta/2 * (Kzp[top] + Kzp[top-1])
        bA[top] =  1 - gamma[top]*dt + beta/2 * (Kzp[top] + Kzp[top-1]) + wsdtdz # okay adding this here 
        dA[top] = Ap[top]

    A = TDMA(aA, bA, cA, dA, self.N)   
    return A


def advance_scalar(self, Kzp, Cp, heat_flux):
    # Initialize tridiagonal arrays for C/temperature. dC is the RHS vector
    aC, bC, cC, dC = initialize_abcd(self.N)
    top = self.top
    beta = self.beta 
    dt  = self.dt 

    aC[1:top] = -0.5*beta*(Kzp[1:top] + Kzp[0:top-1])
    bC[1:top] = 1 + 0.5*beta*(Kzp[2:top+1] + 2*Kzp[1:top] + Kzp[0:top-1])
    cC[1:top] = -0.5*beta*(Kzp[1:top] + Kzp[2:top+1])
    dC[1:top] = Cp[1:top]

    # Bottom-Boundary: no flux for scalars
    bC[0] = 1 + beta/2*(Kzp[1] + Kzp[0])
    cC[0] = -beta/2*(Kzp[1] + Kzp[0])
    dC[0] =  Cp[0] 

    # Top-Boundary: no flux for scalars
    aC[top] = -0.5*beta*(Kzp[top] + Kzp[top-1])
    bC[top] = 1+0.5*beta*(Kzp[top] + Kzp[top-1])
    dC[top] = Cp[-1] + heat_flux # TEMP test
    
    # print('heat flux= ', heat_flux * self.dz )
    C = TDMA(aC, bC, cC, dC, self.N)
    return C

def advance_Q2(self, Q2p, Lp, Kqp, nu_tp, Up, Kzp, N_BV2p, ustar):
    # Initialize tridiagonal arrays for turbulent kinetic energy
    aQ2, bQ2, cQ2, dQ2 = initialize_abcd(self.N)
    top = self.top
    beta = self.beta 
    dt  = self.dt 

    # Dissipation is a size (N-2) vector used for the non-boundary terms in the Q2 equation 
    diss = (2 * dt *(Q2p[1:top]**0.5))/(B1*Lp[1:top])
    aQ2[1:top] = -0.5*beta*(Kqp[1:top] + Kqp[0:top-1])
    bQ2[1:top] = 1 + 0.5*beta*(Kqp[2:top+1] + 2*Kqp[1:top] + Kqp[0:top-1]) + diss 
    cQ2[1:top] = -0.5*beta*(Kqp[1:top] + Kqp[2:top+1])                                  # buoyancy production term  (should be negative such that this term is adding TKE when density is unstable)
    dQ2[1:top] = Q2p[1:top] + 0.25*beta*nu_tp[1:top]*(Up[2:top+1]-Up[0:top-1])**2 - dt*Kzp[1:top]*(N_BV2p[1:top])

    # Bottom-Boundary Condition 
    Q2bot = B1**(2/3) * ustar**2
    bdryterm = 0.5*beta*Kqp[0]*Q2bot
    dissipation = 2 * dt *((Q2p[0]**0.5)/(B1*Lp[0]))
    bQ2[0] = 1+0.5*beta*(Kqp[1] + Kqp[0]) + dissipation
    cQ2[0] = -0.5*beta*(Kqp[1] + Kqp[0])
    dQ2[0] = Q2p[0] + dt*((ustar**4)/nu_tp[0]) - dt*Kzp[0]*(N_BV2p[0]) + bdryterm

    # Top boundary condition
    dissipation =  2 * dt *((Q2p[top]**0.5)/(B1*Lp[top]))
    aQ2[top] = -0.5*beta*(Kqp[top] + Kqp[top-1])
    bQ2[top] = 1+0.5*beta*(Kqp[top] + 2*Kqp[top] + Kqp[top-1]) + dissipation # sw note --> fixed typo w kqp
    dQ2[top] = Q2p[top] + 0.25*beta*nu_tp[top]*((Up[top] - Up[top-1])**2) -4*dt*Kzp[top]*(N_BV2p[top])
    
    Q2  = TDMA(aQ2, bQ2, cQ2, dQ2, self.N) 

    return Q2

def advance_Q2L(self, Q2p, Q2Lp, Lp, Kqp, nu_tp, Up, Kzp, N_BV2p, ustar):

    # Initialize tridiagonal arrays for Q^2 * L (turbulent kinetic energy times a lengthscale)
    aQ2L, bQ2L, cQ2L, dQ2L = initialize_abcd(self.N)
    top = self.top
    beta = self.beta 
    dt  = self.dt 
    z = self.z
    H = self.H

    diss = 2*dt*((Q2p[1:top]**0.5) / (B1*Lp[1:top]))*(1+E2*(Lp[1:top]/(kappa*abs(-H-z[1:top])))**2 \
                                                    + E3*(Lp[1:top]/(kappa*abs(z[1:top])))**2)

    aQ2L[1:top] = -0.5*beta*(Kqp[1:top] + Kqp[0:top-1])
    bQ2L[1:top] = 1 + 0.5*beta*(Kqp[2:top+1] + 2*Kqp[1:top] + Kqp[0:top-1]) + diss
    cQ2L[1:top] = -0.5*beta*(Kqp[1:top] + Kqp[2:top+1]) 
    dQ2L[1:top] = Q2Lp[1:top] + 0.25*beta*nu_tp[1:top]*E1*Lp[1:top] * (Up[2:top+1]-Up[0:top-1])**2 \
                                - 2*dt*Lp[1:top]*E1*Kzp[1:top]*(N_BV2p[1:top])

    # Bottom boundary Condition
    q2lbot = B1**(2/3) * (ustar**2) * kappa * zb
    bdryterm = 0.5*beta*Kqp[0]*q2lbot
    diss =  2 * dt *(Q2p[0]**0.5)/(B1*Lp[0])*(1+E2*(Lp[0]/(kappa*abs(-H-z[0])))**2 + E3*(Lp[0]/(kappa*abs(z[0])))**2)
    bQ2L[0] = 1+0.5*beta*(Kqp[1] + Kqp[0]) + diss
    cQ2L[0] = -0.5*beta*(Kqp[1] + Kqp[0])
    dQ2L[0] = Q2Lp[0] + dt*((ustar**4)/nu_tp[0])*E1*Lp[0] - dt*Lp[0]*E1*Kzp[0]*(N_BV2p[0]) + bdryterm

    # Top boundary condition
    dissipation =  2 * dt *(Q2p[top]**0.5)/(B1*Lp[top])*(1+E2*(Lp[top]/(kappa*abs(-H-z[top])))**2 \
                                                + E3*(Lp[top]/(kappa*abs(z[top])))**2)
    aQ2L[top] = -0.5*beta*(Kqp[top] + Kqp[top-1])
    bQ2L[top] = 1+0.5*beta*(Kqp[top] + 2*Kqp[top] + Kqp[top-1]) + dissipation # Are we using kq or kqp here?
    dQ2L[top] = Q2Lp[top] + 0.25*beta*nu_tp[top]*E1*Lp[top]*(Up[top]-Up[top-1])**2 - 2*dt*Lp[-1]*E1*Kzp[top]*(N_BV2p[top])
    
    Q2L  = TDMA(aQ2L, bQ2L, cQ2L, dQ2L, self.N)
    return  Q2L 