# encoding: utf-8
"""
integrate.py -- A simple RK4 integrator class

Created by Joe Monaco on 2007-11-15.
Copyright (c) 2007 Columbia University. All rights reserved.
"""

class integrator(object):
    
    """
    Per-timestep interface yielding simple RK4 numerical integration. This 
    is not very fast, but it's simple and robust and gives you much more 
    flexibility for coding logic, updating parameters, etc. (You don't have
    to worry about stiff/non-stiff or adaptive timesteps.)
       
       Example usage:
          # initialization
          y = zeros(N)      # initial state
          tau = tau0        # parameter
          D = zeros((m,n))  # parameter
          
          # integration
          myint = integrator(dxdt, y)
          while myint.t < T:
             myint(tau, D)          # integrate the next timestep: y(t+dt) -> y
             y_bar = y.mean()       #   
             tau = new_tau(myint.t) # update parameter
             update(D)              # update parameter
             print 'y_bar(%f) = %f'%(myint.t, y_bar)
          
          # the derivative function
          def dxdt(y, t0, tau, D): 
             return <etc>
    """ 
    
    def __init__(self, deriv, y0, t0=0.0, dt=0.001):
        if not callable(deriv):
            raise TypeError, 'deriv argument must be callable'

        from numpy import array, empty, dot
        
        self.dydt = deriv
        self.dot = dot
        self.t = float(t0)
        self.dt = float(dt)
        self.y = y0
        self._k = empty((4,) + y0.shape, 'd')
        self._c = (self.dt / 6) * array([1, 2, 2, 1], 'd')

    def __call__(self, *args):
        self._k[0,:] = self.dydt(self.y + 0,                        self.t,               *args)
        self._k[1,:] = self.dydt(self.y + self.dt * self._k[0] / 2, self.t + self.dt / 2, *args)
        self._k[2,:] = self.dydt(self.y + self.dt * self._k[1] / 2, self.t + self.dt / 2, *args)
        self._k[3,:] = self.dydt(self.y + self.dt * self._k[2],     self.t + self.dt,     *args)

        self.y += self.dot(self._c, self._k)
        self.t += self.dt

    def next(self, *args): 
        "integrate the next timestep and return the new value"
        self.__call__(*args)
        return self.y
        