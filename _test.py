

def rk4(x, h, y, f):
     k1 = h * f(x, y)
     k2 = h * f(x + 0.5*h, y + 0.5*k1)
     k3 = h * f(x + 0.5*h, y + 0.5*k2)
     k4 = h * f(x + h, y + k3)
     return x + h, y + (k1 + 2*(k2 + k3) + k4)/6.0
    

import numpy

def damped_spring(t, state):
     pos, vel = state
     stiffness = 1
     damping = 0.05
     return numpy.array([vel, -stiffness*pos - damping*vel])

t = 0
dt = 1.0/40
state = numpy.array([5, 0])
print('%10f %10f' % (t, state[0]))

while t < 100:
    t, state = rk4(t, dt, state, damped_spring)
    print('%10f %10f' % (t, state[0]))


