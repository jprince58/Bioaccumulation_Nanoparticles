# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# imports
import numpy as np 
import matplotlib.pyplot as plt 
import nangs
import torch

device = "cuda" if torch.cuda.is_available() else "cpu"

nangs.__version__, torch.__version__

from nangs import PDE

class N6(PDE):
    def computePDELoss(self, inputs, outputs):
                
        #Define outputs
        O=outputs[:,0]
        s=outputs[:,1]
        
        #compute mobile gradients
        grads = self.computeGrads(O, inputs)
        dOdt=grads[:,0]
        dOdx=grads[:,1]
        grads2=self.computeGrads(dOdx, inputs)
        d2Odx2=grads2[:, 1]
        
        #compute attached gradient
        grads = self.computeGrads(s, inputs)
        dsdt=grads[:,0]
        
        #Define independant variables
        t=inputs[:,0]
        x=inputs[:,1]
        
        #Define parameter values
        mu=inputs[:,2]
        nu=inputs[:,3]
        om=inputs[:,4]
        kap=inputs[:,5]
        eps=inputs[:,6]
        rho=inputs[:,7]
        a=inputs[:,8]
        b=inputs[:,9]
        c=inputs[:,10]
        
        #Could I calcualte intermediate values? Will try if it Ig et it to work
        
        # compute loss
        return {
            'mobile': om*(x**b+eps)*d2Odx2 + om*(b*x**(b-1)-2*c*x**(c-1)*(x**b+eps)/(x**c+rho))*dOdx-(om*((c*b*x**((c-1)*(b-1))+c*(c-1)*x**(c-2)*(x**b+eps))/(x**c+rho)-2*c**2*x**(2*c-2)*(x**b+eps)/(x**c+rho)**2)+mu*(nu*(1-x**a)-s)/(x**b+eps))*O+(mu*kap/(x**b+eps))*s-dOdt ,
            'attached': mu*(O*(nu*(1-x**a)-s)-kap*s)/(x**b+eps)-dsdt
            }
    
pde = N6(inputs=('t', 'x', 'mu', 'nu', 'om', 'kap', 'eps', 'rho', 'a', 'b', 'c'), outputs=('O','s'))

from nangs import RandomSampler

sampler = RandomSampler({
    't': [0., 10.], 
    'x': [0., 1.], 
    'mu': [0., 10.],
    'nu': [0., 10.],
    'om': [0., 10.],
    'kap': [0., 10.],
    'eps': [0., 0.75],
    'rho': [0., 0.75],
    'a': [5, 100.],
    'b': [5., 100.],
    'c': [5., 100.],
}, n_samples=1000)

pde.set_sampler(sampler)

# %% Boundary Conditions

x = np.linspace(0, 1, 30)
p0 = np.sin(2*np.pi*x)
p0.shape

from nangs import Dirichlet

n_samples=100

# %%Initial Condition
initial = Dirichlet(
    RandomSampler({
        't': 0.,
        'x': [0., 1.],
        'mu': [0., 10.],
        'nu': [0., 10.],
        'om': [0., 10.],
        'kap': [0., 10.],
        'eps': [0., 0.75],
        'rho': [0., 0.75],
        'a': [5, 100.],
        'b': [5., 100.],
        'c': [5., 100.],
        }, device=device, n_samples=n_samples), 
        lambda inputs: {'O': torch.zeros(n_samples), 's': torch.zeros(n_samples) },
    name="initial"
)

pde.add_boco(initial)

# %%Attached Initial Condition

# initial2 = Dirichlet(
#     RandomSampler({
#         't': 0., 
#         'x': [0., 1.], 
#         'mu': [0., 10.],
#         'nu': [0., 10.],
#         'om': [0., 10.],
#         'kap': [0., 10.],
#         'eps': [0., 0.75],
#         'rho': [0., 0.75],
#         'a': [5, 100.],
#         'b': [5., 100.],
#         'c': [5., 100.],
#         }, device=device, n_samples=1000), 
#         lambda inputs: {'s': 0 },
#     name="initial2"
# )

# pde.add_boco(initial2)

# %%Mobile Supernatant Boundary Condition

BC1 = Dirichlet(
    RandomSampler({
        't': [0., 10.],
        'x': 1., 
        'mu': [0., 10.],
        'nu': [0., 10.],
        'om': [0., 10.],
        'kap': [0., 10.],
        'eps': [0., 0.75],
        'rho': [0., 0.75],
        'a': [5, 100.],
        'b': [5., 100.],
        'c': [5., 100.],
        }, device=device, n_samples=n_samples), 
        lambda inputs: {
            'O': torch.ones(n_samples, device=inputs['x'].device)*(1.+inputs['rho']), 
            's': torch.zeros(n_samples) 
            },
    name="supernatant"
)

pde.add_boco(BC1)

# %%Mobile Substratum Boundary Condition

from nangs import Neumann

class MyNeumann(Neumann):
    def computeBocoLoss(self, inputs, outputs):
        dOdx = self.computeGrads(outputs, inputs)[:, 1]
        return {'gradO': dOdx}

BC2 = MyNeumann(
    RandomSampler({
        't': [0., 10.],
        'x': 0., 
        'mu': [0., 10.],
        'nu': [0., 10.],
        'om': [0., 10.],
        'kap': [0., 10.],
        'eps': [0., 0.75],
        'rho': [0., 0.75],
        'a': [5, 100.],
        'b': [5., 100.],
        'c': [5., 100.],
        }, device=device, n_samples=n_samples), 
    name='substratum'
)

pde.add_boco(BC2)


# %%Define MLP architecture
from nangs import MLP

mlp = MLP(inputs=len(pde.inputs), outputs=len(pde.outputs), layers=5, neurons=200)


# %% Train Model
N_STEPS = 10000

optimizer = torch.optim.Adam(mlp.parameters(), lr=3e-4)
scheduler = torch.optim.lr_scheduler.OneCycleLR(optimizer, max_lr=0.01, pct_start=0.1, total_steps=N_STEPS)
loss_fn = torch.nn.MSELoss()

pde.compile(mlp.to(device), optimizer, scheduler, loss_fn=loss_fn)

hist = pde.solve(N_STEPS)

# %% Plot results
import pandas as pd 

df = pd.DataFrame(hist)
fig = plt.figure(dpi=100)
ax = plt.subplot(1,1,1)
ax.set_yscale('log')
df.plot(ax=ax, grid=True)

t = np.linspace(0,1,50)
x = np.linspace(0,1,100)
mu=1
nu=1
om=1
kap=1
eps=0.2
rho=0.2
a=10
b=10
c=10

grid = np.stack(np.meshgrid(t, x, mu, nu, om, kap, eps, rho, a, b, c), -1).reshape(-1, 2)
X = torch.from_numpy(grid).float().to(device)

p = pde.eval(X).cpu().view(len(t),len(x),1,1,1,1,1,1,1,1,1).numpy()



