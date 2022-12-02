import matplotlib.pyplot as plt
import numpy as np
from numpy import sin as S 
from numpy import cos as C 
from numpy import tan as T
from numpy.core.shape_base import block 
from scipy.linalg import block_diag
import time

class controller:
    def __init__(self):
        self.x, self.y, self.z, self.xdot, self.ydot, self.zdot, self.theta, self.phi, self.psi, self.thetadot, self.psidot, self.phidot = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 
        self.xd, self.yd, self.zd, self.psid = 0.0,0.0,0.0,0.0
        self.g = 9.81
        self.m = 2.0
        self.Ix, self.Iy = 1.2416, 1.2416
        self.Iz = 2.4832
        self.It = np.diag([self.Ix,self.Iy,self.Iz])
        self.F1, self.F2, self.F3, self.F4 = self.m*self.g/4,self.m*self.g/4,self.m*self.g/4,self.m*self.g/4
        self.F1dot, self.F2dot, self.F3dot, self.F4dot = 0.0,0.0,0.0,0.0
        
        self.d = 0.2
        self.c = 0.01
        
        self.kt = np.diag([0.01, 0.01, 0.01])
        self.kr = np.diag([0.001, 0.001, 0.001])
        self.A1,self.A2,self.A3,self.A4,self.A5,self.A6 = np.diag([10,10]),np.diag([10,10]),np.diag([10,10]),np.diag([10,10]),np.diag([1,1]),np.diag([1,1]),
        self.A7 = np.diag([5,5,5,5])
        self.dt = 0.1
        self.t = 0.0

        self.v1,self.v2,self.v3,self.v4,self.v5,self.v6,self.u = np.array([[0.0],[0.0]]),np.array([[0.0],[0.0]]),np.array([[0.0],[0.0]]),np.array([[0.0],[0.0]]),np.array([[0.0],[0.0]]),np.array([[0.0],[0.0]]),np.array([[0.0],[0.0],[0.0],[0.0]])
        self.fx, self.fy, self.fz = 0.0, 0.0, -self.g
        # self.fphi,self.ftheta,self.fpsi = 0.0, 0.0, 0.0
        

    def update(self):
        self.eta = np.array([[self.phi],[self.theta],[self.psi]])
        self.etadot = np.array([[self.phidot],[self.thetadot],[self.psidot]])
        self.g0 = ((self.F1+self.F2+self.F3+self.F4+10e-6)/self.m)*np.array([[S(self.psi), C(self.psi)],[-C(self.psi), S(self.psi)]])
        # print(self.g0)
        self.g1 = np.array([[1/self.Ix , S(self.phi)*T(self.theta)/self.Iy],[0, C(self.phi)/self.Iy]])
        self.g2 = np.array([[C(self.phi)/(self.Iz*C(self.theta)) , 0],[0, C(self.phi)*C(self.theta)/self.m]])

        self.phivec0 = np.array([[S(self.phi)],[C(self.phi)*S(self.theta)]])
        self.phivec1 = np.array([[self.d*(self.F2-self.F4)],[self.d*(self.F3-self.F1)]])
        self.phivec2 = np.array([[self.c*(self.F1-self.F2+self.F3-self.F4)],[self.F1+self.F2+self.F3+self.F4]])

        self.Rt = np.array([[C(self.phi)*C(self.psi), (S(self.phi)*S(self.theta)*C(self.psi) - C(self.phi)*S(self.psi)), (C(self.phi)*S(self.theta)*C(self.psi) + S(self.phi)*S(self.psi))], [C(self.theta)*S(self.psi), (S(self.phi)*S(self.theta)*S(self.psi) + C(self.phi)*C(self.psi)), (C(self.phi)*S(self.theta)*S(self.psi) - S(self.phi)*C(self.psi))], [-S(self.phi), S(self.phi)*C(self.theta), C(self.phi)*C(self.theta)]])
        self.Rr = np.array([[1,0,-S(self.theta)], [0,C(self.phi),C(self.theta)*S(self.phi)], [0,-S(self.phi),C(self.phi)*C(self.theta)]])
        self.dRrbydphi = np.array([[0,0,0], [0,-S(self.phi),C(self.theta)*C(self.phi)], [0,-C(self.phi),-S(self.phi)*C(self.theta)]])
        self.dRrbydtheta = np.array([[0,0,-C(self.theta)], [0,0,-S(self.theta)*S(self.phi)], [0,0,-C(self.phi)*S(self.theta)]])

        # print(np.shape(self.etadot), "etadot")
        # print(np.shape(self.Rr), "Rr")
        # print(np.shape(self.It), "It")
        # print(np.shape(self.Rr.dot(self.etadot)))
        # print(np.shape(self.It.dot(self.Rr).dot(self.etadot)))
        self.farr1 = np.array([[self.fx],[self.fy],[self.fz]])
        self.farr2 = -np.linalg.inv(np.dot(self.It,self.Rr)).dot((self.It.dot(self.phidot*self.dRrbydphi+self.thetadot*self.dRrbydtheta).dot(self.etadot) - np.cross(self.Rr.dot(self.etadot), self.It.dot(self.Rr).dot(self.etadot),axis=0))) + np.array([[(self.c/self.Iz)*C(self.phi)*T(self.theta)*(self.F1-self.F2+self.F3-self.F4)],[(-self.c/self.Iz)*S(self.phi)*(self.F1-self.F2+self.F3-self.F4)],[(self.d/self.Iy)*(S(self.phi)/C(self.theta))*(self.F3-self.F1)]])
        # print(self.farr2)
        # print((self.It.dot(self.phidot*self.dRrbydphi+self.thetadot*self.dRrbydtheta).dot(self.etadot)))# - 
        # print(np.cross(self.Rr.dot(self.etadot), self.It.dot(self.Rr).dot(self.etadot)))
        # print(np.array([[(self.c/self.Iz)*C(self.phi)*T(self.theta)*(self.F1-self.F2+self.F3-self.F4)],[(-self.c/self.Iz)*S(self.phi)*(self.F1-self.F2+self.F3-self.F4)],[(self.d/self.Iy)*(S(self.phi)/C(self.theta))*(self.F3-self.F1)]]))
        # self.farr2 = np.array([[self.fphi],[self.ftheta],[self.fpsi]])
        self.fphi = self.farr2[0]
        self.ftheta = self.farr2[1]
        self.fpsi = self.farr2[2]
        # print(self.ftheta)
        # print(self.fphi)
        # print(self.farr2)

        self.f0 = np.vstack((self.fx,self.fy))
        self.f1 = np.vstack((self.fphi,self.ftheta))
        self.f2 = np.vstack((self.fpsi,self.fz))
        # print(self.f1)
        self.J0 = np.array([[C(self.phi),0],[-S(self.phi)*S(self.theta), C(self.phi)*C(self.theta)]])
        self.J1 = np.array([[0,self.d,0,-self.d],[-self.d,0,self.d,0]])
        self.J2 = np.array([[self.c,-self.c,self.c,-self.c],[1,1,1,1]])

        self.jvec = np.vstack((self.J1,self.J2))
        # print(np.shape(self.jvec))
        # print(np.shape(self.J1), np.shape(self.J2))

        self.gblock = block_diag(self.g1,self.g2)
        # print(self.gblock)

        self.u = np.array([[self.F1dot], [self.F2dot], [self.F3dot], [self.F4dot]])
        # print(self.u)
        # print(np.shape(self.u))

        self.x1 = np.array([[self.x],[self.y]])
        # print(self.x1)
        self.x2 = np.array([[self.xdot], [self.ydot]])
        self.x3 = np.array([[self.phi], [self.theta]])
        self.x4 = np.array([[self.phidot], [self.thetadot]])
        self.x5 = np.array([[self.psi], [self.z]])
        self.x6 = np.array([[self.psidot], [self.zdot]])
        self.x7 = np.array([[self.F1], [self.F2], [self.F3], [self.F4]])
        # print(np.shape(self.x7))

        self.x1dot = self.x2
        self.x2dot = self.f0 + self.g0.dot(self.phivec0)
        self.x3dot = self.x4
        self.x4dot = self.f1 + self.g1.dot(self.phivec1)
        self.x5dot = self.x6
        self.x6dot = self.f2 + self.g2.dot(self.phivec2)
        self.x7dot = self.u
        # print((self.x7*self.dt + self.x7dot)


    def calc_input(self):
        # print(self.x1)
        self.v1 = self.A1@((self.x1d-self.x1)) + self.x1ddot
        # print(self.v1)
        v1dot = self.derivative(self.v1,self.v1o)
        # print(v1dot)
        self.v2 = np.linalg.inv(self.g0)@(((self.x1d - self.x1)) + self.A2@((self.v1 - self.x2)) + v1dot - self.f0 )
        v2dot = self.derivative(self.v2,self.v2o)
        print(np.linalg.inv(self.J0))
        self.v3 = np.linalg.inv(self.J0)@(np.transpose(self.g0)@((self.v1 - self.x2)) + self.A3@(self.v2 - self.phivec0) + v2dot )
        v3dot = self.derivative(self.v3,self.v3o)
        self.v4 = np.linalg.inv(self.g1)@( np.transpose(self.J0)@(self.v2 - self.phivec0) + self.A4@(self.v3 - self.x4) + v3dot - self.f1 )
        v4dot = self.derivative(self.v4,self.v4o)
        self.v5 = self.A5@(self.x5d - self.x5) + self.x5ddot
        v5dot = self.derivative(self.v5,self.v5o)
        self.varr1 = np.vstack((self.v3 - self.x4, self.v5 - self.x6))
        self.v6 = np.linalg.inv(self.g2)@( (self.x5d - self.x5) + self.A6@(self.v5 - self.x6) + v5dot - self.f2 )
        # print(self.v6)
        v6dot = self.derivative(self.v6,self.v6o)
        self.varr2 = np.vstack((self.v4 - self.phivec1,self.v6 - self.phivec2))
        self.u = np.linalg.inv(self.jvec).dot(((np.transpose(self.gblock).dot(self.varr1)) + np.vstack((v4dot,v6dot)) + self.A7.dot(self.varr2)))
        # print(self.u)
        # self.u = np.array([[0.0],[0.0],[0.0],[0.0]])
        self.x7dot = self.u

    def integrate(self,q1,q2):
        return q2*self.dt + q1

    def derivative(self,q1,q2):
        #print("derivaive ", (q1-q2)/self.dt)
        return (q1-q2)/self.dt

    
    def compute(self):
        T = np.linspace(0,20,200)
        self.update()
        
        plt.figure(figsize=(10,5))
        #ax = plt.axes(projection ='3d')

        for i in range(1,200):
            self.t = T[i]
            self.v1o,self.v2o,self.v3o,self.v4o,self.v5o,self.v6o,self.uo = self.v1,self.v2,self.v3,self.v4,self.v5,self.v6,self.u
            #self.update()
            self.x5d = np.array([[0], [5*S(self.t)]])
            self.x5ddot = np.array([[0], [5*C(self.t)]])
            self.x1ddot, self.x1d = np.array([[0], [0]]), np.array([[0], [0]])
            # self.x1ddot, self.x1d = np.array([[0], [0]]), np.array([[0], [0]])
            self.calc_input()

            self.F1, self.F2, self.F3, self.F4 = self.integrate(self.x7,self.x7dot)[:,0]
            # print(self.x7dot,"x7dot")
            # print(self.F1, self.F2, self.F3, self.F4)
            # self.update()
            self.psidot, self.zdot = self.integrate(self.x6,self.x6dot)[:,0]
            # print(self.x6dot,"x6dot")
            # self.update()
            self.psi, self.z = self.integrate(self.x5,self.x5dot)[:,0]
            # print(self.z, "z")
            # self.update()
            self.phidot, self.thetadot = self.integrate(self.x4,self.x4dot)[:,0]
            # print(self.x4, self.x4dot)
            # self.update()
            self.phi, self.theta = self.integrate(self.x3,self.x3dot)[:,0]
            # self.update()
            self.xdot, self.ydot = self.integrate(self.x2,self.x2dot)[:,0]
            # self.update()
            self.x,self.y = self.integrate(self.x1,self.x1dot)[:,0]
            self.update()


            #ax.plot(self.x, self.y, self.z, c='lightblue',marker='o')
            #ax.plot(0, 0, 5*S(self.t), c='red',marker='o')
            plt.plot(self.t,self.z,c='lightblue',marker='o')
            plt.plot(self.t,5*S(self.t),c='red',marker='o')

            #plt.show(block=False)
            plt.pause(0.01)
            # time.sleep(0.1)
        plt.show()

nlc = controller()
nlc.compute()