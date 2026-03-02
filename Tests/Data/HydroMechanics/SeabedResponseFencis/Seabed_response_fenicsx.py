# %% [markdown]
# # Seabed response to peridoc standing water waves

# %%
#FEniCS modules
import dolfinx as dolx
from dolfinx import mesh, fem
from dolfinx.io import XDMFFile
import basix as basix
from dolfinx.fem.petsc import NonlinearProblem
import ufl as ufl

#Python modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpi4py import MPI
import pyvista as pv
import meshio
import os

# %%
dolx.log.set_log_level(dolx.log.LogLevel.INFO)

# %%
print(f"DOLFINx version: {dolx.__version__} based on GIT commit: {dolx.git_commit_hash} of https://github.com/FEniCS/dolfinx/")

# %%
#create a result folder, where all the output is saved
def result_folder(name):
    folder='./'+name

    try: #Checks whether a result folder already exists
        os.mkdir(folder)

    except:
        pass

name_result_folder = 'hm_seabed_response_standing_waves_fenicsx'
result_folder(name_result_folder) #creates a folder for all the results

# %% [markdown]
# ### Analytical solution 2D

# %% [markdown]
# The following 2D analytical solution was written down bei [1]. It was implemented in Python bei Linda Günther [2]. The sign convention follows the sign convention in soil mechanics - pressure is positive and tension negative.

# %%
def compute_pressure_and_stresses(t,x,z):
    
    n=0.4
    ny=0.0              # E = 3K(1-2ny) = 2G(1+ny)
    E = 1e7             # [Pa]
    G=E/(2*(1+ny))      # 100e3    # [Pa]
    K=2/3*G             # [Pa](with ny=0)
   
    Cf=0.001/K          # in the book: Cf = 0.001/K #compressibility of the fluid 
    Cs=0                # compressibility of the solid particles
    Cm=1/K              # compressibility of the porous media
    my=1.3e-3           # [Pa*s] #fluid viscosity of water(Pa s) at T=10°C
    kappa=1.33e-11      # [m²]    (medium sand, kf=10e-4 m/s)   
    gamma_w = 9.81e3    # [Pa/m]
    f_2D = 0.1          # frequency [1/s]
    T_2D = 1/f_2D       # period time [s]
    L_2D = 100          # wave length [m]
    omega=2*np.pi*f_2D  # angular frequency of the wave
    lam=omega*T_2D/L_2D
    
   
    k = kappa*gamma_w/my      #hydraulic conductivity [m/s]    # Gl. (1.33)
    alpha = 1-Cs/Cm                                            # Gl. (4.15)
    S = n*Cf + (alpha-n)*Cs                                    # Gl. (1.28)
    theta = S*G/alpha**2                                       # Gl. (4.13)
    m = 1/(1-2*ny) # = K+1/3*G/G                               # Gl. (4.5)
    cv = k*G*(1+m) / (alpha**2*(1+theta+m*theta)*gamma_w)      # Gl. (4.12)
    xi_2 = complex(lam**2, (omega/cv))                         # Gl. (4.19)
   
    B1 = (1+m)*(xi_2-lam**2)-2*lam*(np.sqrt(xi_2)-lam)
    B2 = 2*m*theta*lam*np.sqrt(xi_2)+theta*((1+m)*(xi_2-lam**2)-2*lam*(np.sqrt(xi_2)-lam))
    B3 = 2*m*theta*lam
    D = 2*lam*(2*lam*(np.sqrt(xi_2)-lam)-(1+m)*(1+m*theta)*(xi_2-lam**2))
    p_rel = np.real((-2*lam*B1*np.exp(-lam*z) - (1+m)*(xi_2-lam**2)*B3*np.exp(-np.sqrt(xi_2)*z))/D * np.exp((omega*t-np.pi*0.5)*1j)*np.cos(lam*x))
    sig_xx_rel = np.real(((-2*(m-1)*lam*theta + 2*lam*(1+m*theta)*lam*z)*B1*np.exp(-lam*z) - 2*lam*B2*np.exp(-lam*z) + ((m-1)*(xi_2-lam**2) - 2*lam**2)*B3*np.exp(-np.sqrt(xi_2)*z))/D * np.exp((omega*t-np.pi*0.5)*1j)*np.cos(lam*x))
    sig_zz_rel = np.real(((-2*(m+1)*lam*theta - 2*lam*(1+m*theta)*lam*z)*B1*np.exp(-lam*z) + 2*lam*B2*np.exp(-lam*z) + ((m-1)*(xi_2-lam**2) + 2*xi_2)*B3*np.exp(-np.sqrt(xi_2)*z))/D * np.exp((omega*t-np.pi*0.5)*1j)*np.cos(lam*x))
    sig_xz_rel = np.real(((-2*lam*(1+m*theta)*lam*z-2*lam*theta)*B1*np.exp(-lam*z) + 2*lam*B2*np.exp(-lam*z) + 2*np.sqrt(xi_2)*lam*B3*np.exp(-np.sqrt(xi_2)*z))/D * np.exp((omega*t-np.pi*0.5)*1j)*np.sin(lam*x))
    return p_rel, sig_xx_rel, sig_zz_rel, sig_xz_rel

# %%
L_2D = 100
y_2D = np.linspace(0,100,1000)
y_rel = y_2D/L_2D
colors = {0:"orangered", 2:"gold", 4:"blueviolet", 6:"forestgreen", 8:"darkorange", 10:"royalblue"}

fig, ax = plt.subplots(ncols=2, figsize=(15,7))
ax[0].spines['right'].set_color('none') #https://www.grund-wissen.de/informatik/python/scipy/matplotlib.html
ax[0].spines['top'].set_color('none') #https://www.grund-wissen.de/informatik/python/scipy/matplotlib.html
ax[1].spines['right'].set_color('none') #https://www.grund-wissen.de/informatik/python/scipy/matplotlib.html
ax[1].spines['top'].set_color('none') #https://www.grund-wissen.de/informatik/python/scipy/matplotlib.html

for idx in (0,1):
    ax[idx].grid(True)
    ax[idx].set_ylabel("$y$ / $L$")
    ax[idx].set_xlim(-1.1,1.1)

for t in [0,2,4,6,8,10]:
    ax[0].plot(compute_pressure_and_stresses(t,0,y_2D)[0], -y_rel, color=colors[t], label= "t = %.1f s" %t)

t=2.5
ax[0].set_xlabel("$p$ / $\\tilde{p}$")
ax[0].legend()
ax[1].plot(compute_pressure_and_stresses(t,0,y_2D)[1], -y_rel, color = colors[6], label = r"$\sigma'_{xx}$ / $(\alpha\tilde{p})$")
#ax[1].plot(compute_pressure_and_stresses(t,0,y)[1]+compute_pressure_and_stresses(t,0,y)[0], -y_rel, linestyle = "--", color = colors[3], label = "$\\sigma_{xx}$/$\\alpha\\tilde{p}$")  # Total horizontal stress
ax[1].plot(compute_pressure_and_stresses(t,0,y_2D)[2], -y_rel, color = colors[2], label = r"$\sigma'_{yy}$ / $(\alpha\tilde{p})$")
#ax[1].plot(compute_pressure_and_stresses(t,0,y)[2]+compute_pressure_and_stresses(t,0,y)[0], -y_rel, linestyle = "--", color = colors[1], label = "$\\sigma_{yy}$/$\\alpha\\tilde{p}$")  # Total vertical stress
ax[1].plot(compute_pressure_and_stresses(t,0,y_2D)[3], -y_rel, color = colors[4], label = r"$\sigma'_{xy}$ / $(\alpha\tilde{p})$")
ax[1].set_xlabel(r"$\sigma'$ / $(\alpha\tilde{p})$")
ax[1].legend();
plt.savefig('./'+name_result_folder+'/distribution_over_depth_over_time_analytical_2D.pdf', bbox_inches='tight', format = 'pdf', dpi=600) 

# %%
t=np.linspace(0,20,2000)
colors = {1:"gold", 2:"blueviolet", 3:"forestgreen", 4:"royalblue"}

fig, ax = plt.subplots()
ax.spines['right'].set_color('none') #https://www.grund-wissen.de/informatik/python/scipy/matplotlib.html
ax.spines['top'].set_color('none') #https://www.grund-wissen.de/informatik/python/scipy/matplotlib.html

#fig, ax = plt.subplots(ncols=2, figsize=(15,7))
#for idx in (0,1):
#    ax[idx].grid(True)
#    ax[idx].set_xlabel("$t$ / s")
ax.set_xlabel("$t$ / s")
for y in (np.linspace(0,100,6)):
    ax.plot( t, compute_pressure_and_stresses(t,0,y)[0], color = colors[4])
    ax.set_ylabel("$p$ / $\\tilde{p}$")
   # ax[1].plot(t, compute_pressure_and_stresses(t,0,y)[1], color = colors[3], label = r"$\sigma'_{xx}/(\alpha\tilde{p})$")
   # ax[1].plot(t, compute_pressure_and_stresses(t,0,y)[2], color = colors[1], label = r"$\sigma'_{yy}/(\alpha\tilde{p})$")
   # ax[1].plot(t, compute_pressure_and_stresses(t,0,y)[3], color = colors[2], label = r"$\sigma'_{xy}/(\alpha\tilde{p})$")
   # if y == 0:
   #     ax[1].legend(loc="upper right")

#ax[1].set_ylabel("$\sigma$'/$\\alpha\\tilde{p}$")

    
ax.set_title("Pore water pressure over time")
plt.savefig('./'+name_result_folder+'/pore_water_pressure_over_time_analytical_2D.pdf', bbox_inches='tight', format = 'pdf', dpi=600) 
#ax[1].set_title("Effective stresses over time");

# %%
x, y = np.meshgrid(np.linspace(0,200,1000),np.linspace(0,100,1000))
t = 2.5

fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(15,7))
l1=ax[0][0].contourf(x,-y, compute_pressure_and_stresses(t,x,y)[0], 15)
l2=ax[0][1].contourf(x,-y, compute_pressure_and_stresses(t,x,y)[1], 15)
l3=ax[1][1].contourf(x,-y, compute_pressure_and_stresses(t,x,y)[2], 15)
l4=ax[1][0].contourf(x,-y, compute_pressure_and_stresses(t,x,y)[3], 15)
fig.colorbar(l1,ax=ax[0][0])
fig.colorbar(l2,ax=ax[0][1])
fig.colorbar(l3,ax=ax[1][1])
fig.colorbar(l4,ax=ax[1][0])
for i in (0,1):
    for j in (0,1):
        ax[i][j].set_aspect('equal')
        ax[i][j].set_xlabel('$x$ / m')
        ax[i][j].set_ylabel('$y$ / m')
ax[0][0].set_title("$p$ / $\\tilde{p}$")
ax[0][1].set_title("$\\sigma'_{xx}$ / $\\alpha\\tilde{p}$")
ax[1][1].set_title("$\\sigma'_{yy}$ / $\\alpha\\tilde{p}$")
ax[1][0].set_title("$\\sigma'_{xy}$ / $\\alpha\\tilde{p}$")
fig.tight_layout();
plt.savefig('./'+name_result_folder+'/distribution_over_depth_analytical_2D.pdf', bbox_inches='tight', format = 'pdf', dpi=600) 

# %% [markdown]
# ## FEniCsx

# %% [markdown]
# ### Mesh generation and definition of boundary regions

# %%
height = 100 #m #initial height of the seabed
width  = 200 #m
Ny = 25
Nx = 50

domain = mesh.create_rectangle(MPI.COMM_WORLD,[np.array([0,0]), np.array([width, height])], [Nx, Ny], mesh.CellType.triangle)
cell_name = domain.topology.cell_name() 
gdim = domain.geometry.dim

tdim = domain.topology.dim
degree = 1
shape = (tdim,) #vector field of size dim

fdim = domain.topology.dim - 1 #facet dimension
vdim = domain.topology.dim - 2 #vertices dimension

Vue = basix.ufl.element("Lagrange", cell_name, degree=2, shape=(gdim,))  # displacement finite element
Vpe = basix.ufl.element("Lagrange", cell_name, degree=1)                 # pore pressure finite element
mixed = basix.ufl.mixed_element([Vue,Vpe])
V = fem.functionspace(domain, mixed)

n_reference_fa3 = ufl.as_vector((0,1)) #normal vector reference configuration top facet
dimension = '2D'

# %%
pv.set_plot_theme("document")
pv.set_jupyter_backend("static")

topology, cell_types, geometry = dolx.plot.vtk_mesh(domain, domain.geometry.dim)
grid = pv.UnstructuredGrid(topology, cell_types, geometry)

plotter = pv.Plotter()
plotter.add_mesh(grid, show_edges=True)
plotter.view_xy()
if not pv.OFF_SCREEN:
    plotter.show()
else:
    figure = plotter.screenshot("mesh.png")

# %%
#define boundary regions
def bottom(x):
    return np.isclose(x[1], 0) 

def top(x):
    return np.isclose(x[1], height) 

def left(x):
    return np.isclose(x[0], 0) 

def right(x):
    return np.isclose(x[0], width)

bottom_facets   = mesh.locate_entities_boundary(domain, fdim, bottom)
top_facets      = mesh.locate_entities_boundary(domain, fdim, top)
left_facets     = mesh.locate_entities_boundary(domain, fdim, left)
right_facets    = mesh.locate_entities_boundary(domain, fdim, right)

# %%
#Concatenate and sort the arrays based on facet indices. bottom facets marked with 1, right facets with two,...
marked_facets = np.hstack([bottom_facets, right_facets, top_facets,left_facets])
marked_values = np.hstack([np.full_like(bottom_facets, 1), np.full_like(right_facets, 2),np.full_like(top_facets, 3),np.full_like(left_facets, 4)])
sorted_facets = np.argsort(marked_facets)
facet_tag     = mesh.meshtags(domain, fdim, marked_facets[sorted_facets], marked_values[sorted_facets])

# %% [markdown]
# ### Material characteristics

# %%
def lame(E, nu):  # [1]
    mu = E / (2 * (1 + nu))
    lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))  # plane strain
    return mu, lmbda

# %%
E           = fem.Constant(domain,1e7)      #N/m² Youngs Modulus
nu          = fem.Constant(domain,0.0)      #PoissonRatio
mu, lmbda   = lame(E,nu)

rho_sr      = fem.Constant(domain,2650.0)   #kg/m³  #density of the soil particles
rho_lr      = fem.Constant(domain,1000.0)   #kg/m³  #density of water (the fluid)

n_f_ini     = fem.Constant(domain,0.4)      #initial porosity -volume fraction of the fluid

k_s         = ufl.as_tensor([[1.33e-11,0],[0,1.33e-11]]) #m² permeability
mu_fr       = fem.Constant(domain,1.3e-3)   #fluid viscosity of water(Pa s) at T=10°C

#Parameter for the fluid
beta_m      = fem.Constant(domain,4.81e-10) #1/Pa fluid compressibility
p_ref       = fem.Constant(domain,101325.0) #Pa #reference pressure - in this case standard pressure
rho_fr_ref  = fem.Constant(domain,1000.0)

b  = fem.Constant(domain,(0.0, 0.0))
model = "seabed"

# %% [markdown]
# ### Time period

# %%
# defines the time steps dt and the overall time limit
t_      = 0.0
t_end   = 10
Nsteps  = 100
dt_     = t_end/Nsteps
dt      = fem.Constant(domain,dt_) 

# %% [markdown]
# ### Functions

# %%
U       = fem.Function(V)   #the whole solution containing displacement and porewater pressure
(u, p)  = ufl.split(U)      #splitting the solution into displacement u and porewater pressure p

# %%
v           = ufl.TestFunction(V)   #the whole test function
(vu, vp)    = ufl.split(v)          #splitting the whole test function into two parts

dU          = ufl.TrialFunction(V)  #the whole trial function
(du, dp)    = ufl.split(dU)         #splitting the whole trial function into two parts

# %%
ram_steps = 5              #load ramp in the first 5 steps

# %% [markdown]
# ### Dirichlet Boundary conditions

# %%
amplitude_wave  = 0.1e5 #Pa
f               = 0.1   #1/s
w               = 2*np.pi*f  #1/s
L               = 100   #m

A_t             = fem.Constant(domain, 0.0)

# %%
#Dirichlet Boundary Conditions

#hydraulic boundary conditions concerning the pore water pressure
V_p, mapping = V.sub(1).collapse()  
p_wave = fem.Function(V_p)

top_dofs_p = fem.locate_dofs_topological((V.sub(1),V_p),facet_tag.dim, facet_tag.find(3))

#mechanical boundary condition concerning the displacement
V_u, mapping = V.sub(0).collapse()

V_ux, mapping = V.sub(0).sub(0).collapse() 
u_x = fem.Function(V_ux)

V_uy, mapping = V.sub(0).sub(1).collapse() 
u_y = fem.Function(V_uy)

bottom_dofs_uy = fem.locate_dofs_topological((V.sub(0).sub(1),V_uy), facet_tag.dim, facet_tag.find(1))
right_dofs_ux  = fem.locate_dofs_topological((V.sub(0).sub(0),V_ux), facet_tag.dim, facet_tag.find(2))
left_dofs_ux   = fem.locate_dofs_topological((V.sub(0).sub(0),V_ux), facet_tag.dim, facet_tag.find(4))

bc_h1 = fem.dirichletbc(p_wave, top_dofs_p,V.sub(1))           #hydraulic wave boundary
bc_m1 = fem.dirichletbc(u_y, bottom_dofs_uy, V.sub(0).sub(1))  #no movement in y direction of the bottom facet
bc_m2 = fem.dirichletbc(u_x, right_dofs_ux, V.sub(0).sub(0))   #no movement in x direction of the right facet
bc_m3 = fem.dirichletbc(u_x, left_dofs_ux, V.sub(0).sub(0))    #no movement in x direction of the left facet

# %% [markdown]
# ### Neumann Boundary conditions

# %%
def kinematics(u):
    d = len(u)
    I = ufl.Identity(d)             # Identity tensor
    F = I + ufl.grad(u)             # Deformation gradient
    J = ufl.det(F)                  # Jacobian of the deformation gradient 
    cof_F = J*ufl.inv(F.T)
    b_s = ufl.dot(F,F.T)
    b_inv = ufl.inv(b_s)
    C = ufl.dot(F.T,F)  		    # Right Cauchy Green tensor
    C_inv = ufl.inv(C)
    return d, I, F, J, cof_F, b_s,b_inv,C, C_inv

# %%
def neumann_condition(dimension,nc,n_reference_fa,p_fa):
    n_reference = n_reference_fa #normal vector in the reference configuration
    if nc == 'constant': # load stays the same and keep its oritentation
        T_fa = 1*p_fa*n_reference
    else:
        print('nc not defnined')

    if dimension == '2D':
        T = ufl.as_vector((T_fa[0],T_fa[1]))
    else:
        T = T_fa
        
    return T

# %%
#type and value of Neumann Boundary Condition
nc      = 'constant'
x_ufl   = ufl.SpatialCoordinate(domain)
p_fa3   = -A_t*ufl.cos(((2*ufl.pi)/L)*x_ufl[0])

# %%
Q  = fem.Constant(domain,0.0) # Fluid flow on the boundary

# %% [markdown]
# ### Initial conditions

# %%
#poromechanics initial conditions
Uold = fem.Function(V, name="Initial")
uold, pold = ufl.split(Uold)

u_derivative = (u - uold) / dt  
p_derivative = (p - pold) / dt

# %% [markdown]
# ### Variational formulation

# %%
def porosity_development(u,n_f_ini):
    d, I, F_s, J_s,cof_F_s,b_s,b_inv, C_s, C_s_inv = kinematics(u)
    n_f = 1-((1-n_f_ini)/J_s)   #volume fraction fluid 
    e = n_f/(1-n_f)
    n_s = 1-n_f                 #saturation condition 
    return n_f, e, n_s

def density_development(u,n_f_ini,rho_sr,rho_fr):
    n_f, e, n_s = porosity_development(u,n_f_ini)
    rho_s = n_s*rho_sr
    rho_f = n_f*rho_fr
    density = rho_s + rho_f
    return density

# %%
def Second_piola_kirchhff(u,mu,lmbda):
    d, I, F, J,cof_F,b_s,b_inv, C, C_inv = kinematics(u)
    S = mu*(I-C_inv)+lmbda*ufl.ln(J)*C_inv #derived from the energy_function for compressible Neo-Hooke
    return S

def cauchy(u,mu, lmbda):
    d, I, F, J,cof_F, b_s,b_inv,C, C_inv = kinematics(u)
    S = Second_piola_kirchhff(u,mu,lmbda)
    sigma = ufl.inv(J)*ufl.dot(ufl.dot(F,S),F.T)
    return sigma

# %%
def d_green_lagrange(u,v):
    grad_v = ufl.grad(v)
    grad_u = ufl.grad(u)
    dEgl =(1/2)*(grad_v+grad_v.T+ufl.dot(grad_u.T,grad_v)+ufl.dot(grad_v.T,grad_u))
    return dEgl

# %%
def variational_form_hm_mech(u,p,vu,mu,lmbda,n_f_ini,rho_sr,rho_fr_ref,b,p_fa3,nc,n_reference_fa3,dimension,dx,ds):
    d, I, F_s, J_s,cof_F_s,b_s,b_inv, C_s, C_s_inv = kinematics(u)
    rho_fr = rho_fr_ref*ufl.exp(beta_m*(p-p_ref))
    rho = density_development(u, n_f_ini, rho_sr,rho_fr)
    rho_0 = J_s*rho
    S_eff = Second_piola_kirchhff(u,mu,lmbda)
    T_fa3 = neumann_condition(dimension,nc,n_reference_fa3, p_fa3)
    dEgl_term = d_green_lagrange(u,vu)
    vol_body_force = rho_0*b #volumetric body force
    mech_form = ufl.inner(S_eff, dEgl_term)*dx - p*ufl.inner(J_s*C_s_inv, dEgl_term)*dx - ufl.inner(vol_body_force, vu)*dx - ufl.inner(T_fa3,vu)*ds(3)
    return mech_form

def variational_form_hm_hydro(u,u_derivative,p,p_derivative,vp,beta_m,rho_fr_ref,p_ref,k_s,mu_fr,n_f_ini,b,Q,dx,ds):
    d,I, F_s, J_s,cof_F_s, b_s,b_inv,C_s, C_s_inv = kinematics(u)
    n_f,e,n_s = porosity_development(u,n_f_ini) # n_f and n_s as current volume fractions
    rho_fr = rho_fr_ref*ufl.exp(beta_m*(p-p_ref))      #pressure-dependent density
    k = k_s/mu_fr
    K = J_s*ufl.dot(ufl.dot(ufl.inv(F_s),k),ufl.inv(F_s.T))
    grad_p = ufl.grad(p)
    grad_vp = ufl.grad(vp)
    dEgl = d_green_lagrange(u,u_derivative)
    hydro_form = rho_fr*ufl.inner(J_s*C_s_inv, dEgl)*vp*dx +rho_fr*ufl.dot(ufl.dot(grad_vp,K),grad_p)*dx-(rho_fr**2)*ufl.dot(grad_vp,ufl.dot(K,(ufl.dot(b,F_s))))*dx + rho_fr*J_s*n_f*beta_m*p_derivative*vp*dx- rho_fr*ufl.inner(Q,vp)*ds(3)
    return hydro_form

# %%
ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tag)
dx = ufl.Measure("dx",domain=domain, metadata={"quadrature_degree":4}) 
mech_form  = variational_form_hm_mech(u,p,vu,mu,lmbda,n_f_ini,rho_sr,rho_fr_ref,b,p_fa3,nc,n_reference_fa3,dimension,dx,ds)
hydro_form = variational_form_hm_hydro(u,u_derivative,p,p_derivative,vp,beta_m,rho_fr_ref,p_ref,k_s,mu_fr,n_f_ini,b, Q,dx,ds)
FF         = mech_form+hydro_form

# %%
petsc_options = {

    "snes_type": "newtonls",                  # Nonlinear Solver
    "snes_linesearch_type": "none",
    "snes_rtol": 1e-9,
    "snes_monitor": None,
   
    "pc_type": "lu",                          # Preconditioner type
    "pc_factor_mat_solver_type": "mumps",     # MUMPS direkt

    "ksp_type": "preonly",                    # Linear solver type
    "ksp_rtol": 1.0e-18,
    "ksp_atol": 2e-7,
    "ksp_max_it": 2000,
    "ksp_monitor": None
}

# %%
# Compute Jacobian of FF 
JC = ufl.derivative(FF,U,dU)

# %%
problem = NonlinearProblem(
    FF,
    U,
    bcs=[bc_h1,bc_m1,bc_m2,bc_m3],
    petsc_options=petsc_options,
    petsc_options_prefix="seabed_",
    J= JC)

# %%
file_results_u      = XDMFFile(domain.comm, "./"+name_result_folder+"/seabed_"+nc+"_u.xdmf", "w")
file_results_p      = XDMFFile(domain.comm, "./"+name_result_folder+"/seabed_"+nc+"_p.xdmf", "w")
file_results_stress = XDMFFile(domain.comm, "./"+name_result_folder+"/seabed_"+nc+"_stress.xdmf", "w")
file_results_u.write_mesh(domain)
file_results_p.write_mesh(domain)
file_results_stress.write_mesh(domain)

# %% [markdown]
# ### Result preparation

# %%
gdim = domain.geometry.dim
k = V.ufl_element().degree
V_stress = fem.functionspace(domain,("CG",k-1, (gdim,gdim))) # DG1
stress_expr = cauchy(u,mu,lmbda) #symbolic UFL form
expr = fem.Expression(stress_expr, V_stress.element.interpolation_points) 

stress_h = fem.Function(V_stress, name="Cauchy") 
stress_h.interpolate(expr) #real numbers

# %% [markdown]
# ### Solving

# %%
V_u_out = fem.functionspace(domain,basix.ufl.element("Lagrange",cell_name,degree=1,shape=(gdim,)))
V_p_out = fem.functionspace(domain,basix.ufl.element("Lagrange",cell_name,degree=1))

u_out = fem.Function(V_u_out, name="Displacement")
p_out = fem.Function(V_p_out, name="Porewater pressure")

# %%
time = np.linspace(0,t_end,Nsteps+1) #[5]
i = 0
for n, t_ in enumerate(time): 
    print("Steps:",i)
    print("Time:", t_)
    
    ramp_factor = min(n/ram_steps, 1.0)
    A_t.value = ramp_factor * amplitude_wave * np.sin(w*t_)

    p_wave.interpolate(
        lambda x: A_t.value * np.cos(((2*np.pi)/L)*x[0])
    )
    
    problem.solve()
    converged = problem.solver.getConvergedReason()
    num_iter = problem.solver.getIterationNumber()
    assert converged > 0, f"Solver did not converge, got {converged}."
    print(
        f"Solver converged after {num_iter} iterations with converged reason {converged}."
    )

    u_out.interpolate(U.sub(0))
    p_out.interpolate(U.sub(1))

    stress_h.interpolate(expr)

    file_results_u.write_function(u_out, t_)
    file_results_p.write_function(p_out, t_)
    file_results_stress.write_function(stress_h, t_)

    Uold.x.array[:] = U.x.array
    Uold.x.scatter_forward()

    i = i+1

file_results_u.close()
file_results_p.close()
file_results_stress.close()

# %% [markdown]
# ## Post processing

# %%
def xdmf_to_vtu(xdmf_path,step):   
    with meshio.xdmf.TimeSeriesReader(xdmf_path+".xdmf") as reader:
    
        points, cells = reader.read_points_cells()
    
        if points.shape[1] == 2:
            points = np.column_stack([points, np.zeros(points.shape[0])])

    
        for k in range(0,reader.num_steps,step):

            t, point_data, cell_data = reader.read_data(k)

            mesh = meshio.Mesh(
                points=points,
                cells=cells,
                point_data=point_data,
                cell_data=cell_data
            )

            outname = f"{xdmf_path}_step_{k:04d}_t_{t:.3f}.vtu"
            meshio.write(outname, mesh)

# %%
xdmf_to_vtu('./'+name_result_folder+'/seabed_constant_p',20)
xdmf_to_vtu('./'+name_result_folder+'/seabed_constant_stress',20)

# %%
reader = pv.get_reader('./'+name_result_folder+'/seabed_constant_p_step_0000_t_0.000.vtu')
mesh_t0 = reader.read()

reader = pv.get_reader('./'+name_result_folder+'/seabed_constant_p_step_0020_t_2.000.vtu')
mesh_t2 = reader.read()

reader = pv.get_reader('./'+name_result_folder+'/seabed_constant_p_step_0040_t_4.000.vtu')
mesh_t4 = reader.read()

reader = pv.get_reader('./'+name_result_folder+'/seabed_constant_p_step_0060_t_6.000.vtu')
mesh_t6 = reader.read()

reader = pv.get_reader('./'+name_result_folder+'/seabed_constant_p_step_0080_t_8.000.vtu')
mesh_t8 = reader.read()

reader = pv.get_reader('./'+name_result_folder+'/seabed_constant_p_step_0100_t_10.000.vtu')
mesh_t10 = reader.read()

plotter = pv.Plotter()
sargs = dict(title="p / Pa", height=0.3, position_x=0.2, position_y=0.02)  # noqa: C408
plotter.add_mesh(
    mesh_t4,
    scalars="Porewater pressure",
    show_edges=False,
    show_scalar_bar=True,
    label="p",
    scalar_bar_args=sargs,
)
plotter.show_bounds(ticks="outside", xtitle="x / m", ytitle="y / m")
plotter.add_axes()
plotter.view_xy()
plotter.show()

# %%
reader = pv.get_reader('./'+name_result_folder+'/seabed_constant_stress_step_0000_t_0.000.vtu')
mesh_t0_stress = reader.read()

reader = pv.get_reader('./'+name_result_folder+'/seabed_constant_stress_step_0020_t_2.000.vtu')
mesh_t2_stress = reader.read()

reader = pv.get_reader('./'+name_result_folder+'/seabed_constant_stress_step_0040_t_4.000.vtu')
mesh_t4_stress = reader.read()

reader = pv.get_reader('./'+name_result_folder+'/seabed_constant_stress_step_0060_t_6.000.vtu')
mesh_t6_stress = reader.read()

reader = pv.get_reader('./'+name_result_folder+'/seabed_constant_stress_step_0080_t_8.000.vtu')
mesh_t8_stress = reader.read()

reader = pv.get_reader('./'+name_result_folder+'/seabed_constant_stress_step_0100_t_10.000.vtu')
mesh_t10_stress = reader.read()

# %%
colors = {
    0: "orangered",
    2: "gold",
    4: "blueviolet",
    6: "forestgreen",
    8: "darkorange",
    10: "royalblue",
}

# %%
def nodes_along_line(mesh):
    tolerance = 1e-8
    mask = np.abs(mesh.points[:,0]) < tolerance     #boolean whether node lays at x=0
    line_points = mesh.points[mask]
    order = np.argsort(line_points[:,1])
    line_points = line_points[order]                #node points along x=0
    
    return mask, line_points

def get_pressure_sorted(mesh,mask,line_points):
    pressure = mesh.point_data["Porewater pressure"][mask]
    depth = line_points[:, 1]
    indices_sorted = np.argsort(depth)
    return pressure[indices_sorted]

def get_depth_sorted_nodes(line_points):
    depth = line_points[:, 1]
    indices_sorted = np.argsort(depth)
    return depth[indices_sorted]

def get_stresses_sorted(mesh, mask, line_points):
    sigma = mesh.point_data["Cauchy"][mask]
    depth = line_points[:,1]
    indices_sorted = np.argsort(depth)
    sigma_xx = -sigma[indices_sorted,0]
    sigma_yy = -sigma[indices_sorted,3]
    sigma_xy = -sigma[indices_sorted,1]
    return sigma_xx,sigma_yy,sigma_xy

# %%
x = 0
y = np.linspace(0, 100, 1000)
y_rel = y / 100

fig, ax = plt.subplots(figsize=(9, 6))
ax.spines['right'].set_color('none') #https://www.grund-wissen.de/informatik/python/scipy/matplotlib.html
ax.spines['top'].set_color('none') #https://www.grund-wissen.de/informatik/python/scipy/matplotlib.html

## Plotting analytical solution
for t in [2, 4, 6, 8, 10]:
    if t !=10:
        ax.plot(
            compute_pressure_and_stresses(t, x, y)[0],
            -y_rel,
            color=colors[t],
            alpha = 0.6,
        )
    else:
        ANALYTICAL, = ax.plot(
            compute_pressure_and_stresses(t, x, y)[0],
            -y_rel,
            color=colors[t],
            alpha = 0.6,
            label=f"analytical",
        )

## Plotting numerical solution
i_steps=[2,4,6,8,10]
meshes = [mesh_t2,mesh_t4,mesh_t6,mesh_t8,mesh_t10]

for step, mesh in zip(i_steps, meshes):
    mask, line_points = nodes_along_line(mesh)
    pressure = get_pressure_sorted(mesh,mask,line_points)
    depth = get_depth_sorted_nodes(line_points)
    if step != 10:
        ax.plot(
            pressure/0.1e5,
            (depth/100)-1,
            'o',
            color=colors[step]
            )
        ax.set_xlabel("$p$ / $\\tilde{p}$")
    else:
        ax.plot(
            pressure/0.1e5,
            (depth/100)-1,
            'o',
            color=colors[step],
            label=f"numerical - FEniCSx"
            )
        ax.set_xlabel("$p$ / $\\tilde{p}$")


# Add the first legend
first_patch = mpatches.Patch(color=colors[2], label='$t = 2.0$ s')
second_patch = mpatches.Patch(color=colors[4], label='$t = 4.0$ s')
third_patch = mpatches.Patch(color=colors[6], label='$t = 6.0$ s')
forth_patch = mpatches.Patch(color=colors[8], label='$t = 8.0$ s')
fifth_patch = mpatches.Patch(color=colors[10], label='$t = 10.0$ s')

first_legend = ax.legend(handles=[first_patch, second_patch, third_patch,forth_patch,fifth_patch], loc='lower left')

# Manually add the first legend back to the plot
ax.add_artist(first_legend),

ax.set_ylabel("$y$ / $L$")
ax.set_xlim(-1.1, 1.1)
ax.legend(bbox_to_anchor=(.65, 0.15))
plt.savefig('./'+name_result_folder+'/comparison_analytical_fenicsx_p.pdf')

# %%
x = 0
y = np.linspace(0, 100, 1000)
y_rel = y / 100

fig, ax = plt.subplots(figsize=(9, 6))
ax.spines['right'].set_color('none') #https://www.grund-wissen.de/informatik/python/scipy/matplotlib.html
ax.spines['top'].set_color('none') #https://www.grund-wissen.de/informatik/python/scipy/matplotlib.html

## Plotting analytical solution
for t in [2, 4, 6, 8, 10]:
    if t !=10:
        ax.plot(
            compute_pressure_and_stresses(t, x, y)[2],
            -y_rel,
            color=colors[t],
            alpha = 0.6,
        )
    else:
        ANALYTICAL, = ax.plot(
            compute_pressure_and_stresses(t, x, y)[2],
            -y_rel,
            color=colors[t],
            alpha = 0.6,
            label=f"analytical",
        )


## Plotting numerical solution
i_steps=[2,4,6,8,10]
meshes_stress = [mesh_t2_stress,mesh_t4_stress,mesh_t6_stress,mesh_t8_stress,mesh_t10_stress]

for step, mesh in zip(i_steps, meshes_stress):
    mask,line_points = nodes_along_line(mesh)
    sigma_xx, sigma_yy, sigma_xy = get_stresses_sorted(mesh,mask,line_points)
    depth =  get_depth_sorted_nodes(line_points)

    if step != 10:
        ax.plot(
            sigma_yy/0.1e5,
            (depth / 100)-1,
            'o',
            color=colors[step]
            )
        ax.set_xlabel("$\\sigma'_{yy,\\text{FEniCSx}}$ / $\\alpha\\tilde{p}$")
    else:
        ax.plot(
            sigma_yy/0.1e5,
            (depth / 100)-1,
            'o',
            color=colors[step],
            label=f"numerical - FEniCSx"
            )
        ax.set_xlabel("$\\sigma'_{yy,\\text{FEniCSx}}$ / $\\alpha\\tilde{p}$")


# Add the first legend
first_patch = mpatches.Patch(color=colors[2], label='$t = 2.0$ s')
second_patch = mpatches.Patch(color=colors[4], label='$t = 4.0$ s')
third_patch = mpatches.Patch(color=colors[6], label='$t = 6.0$ s')
forth_patch = mpatches.Patch(color=colors[8], label='$t = 8.0$ s')
fifth_patch = mpatches.Patch(color=colors[10], label='$t = 10.0$ s')

first_legend = ax.legend(handles=[first_patch, second_patch, third_patch,forth_patch,fifth_patch], loc='lower left')

# Manually add the first legend back to the plot
ax.add_artist(first_legend),

ax.set_ylabel("$y$ / $L$")
ax.legend(bbox_to_anchor=(.65, 0.15))
plt.savefig('./'+name_result_folder+'/comparison_analytical_fenicsx_stress.pdf')

# %%
def compute_abs_and_rel_pressure_error(pressures, depth, t, x):
    num_points = pressures.shape[0]
    f_abs = np.zeros(num_points)
    
    analytical_index = np.arange(num_points)
    numerical_index  = np.arange(num_points-1,-1,-1)

    for ana_idx,num_idx in zip(analytical_index, numerical_index):
        y = depth[ana_idx]
        pressure_ana = compute_pressure_and_stresses(t, x, y)[0]  # returns pressure normalised to the pressure amplitude #calculated from top to bottom
        pressure_num = (pressures[num_idx]/0.1e5)                 # absolute pressure divided by pressure amplitude       #calculated originally from bottom to top - num idx guarantees that parameters from top to bottom are called
        
        f_abs[num_idx] = pressure_num - pressure_ana              #saving with num_idx guarantees a saving from bottom to top (later it is plotted from bottom to top)


    return f_abs

def compute_abs_and_rel_stress_error(sigma_yy,sigma_xx,sigma_xy, depth, t, x):
    num_points = depth.shape[0]
    f_abs_sigma_yy = np.zeros((num_points))
    f_abs_sigma_xx = np.zeros((num_points))
    f_abs_sigma_xy = np.zeros((num_points))
    analytical_index = np.arange(num_points)
    numerical_index  = np.arange(num_points-1,-1,-1)

    for ana_idx,num_idx in zip(analytical_index, numerical_index):
        y = depth[ana_idx]
        sigma_ana_yy = compute_pressure_and_stresses(t, x, y)[2]  # returns stresses normalised to the pressure amplitude   #calculated from top to bottom
        sigma_ana_xx = compute_pressure_and_stresses(t, x, y)[1]  # returns stresses normalised to the pressure amplitude
        sigma_ana_xy = compute_pressure_and_stresses(t, x, y)[3]
        
        sigma_num_yy = (sigma_yy[num_idx]/ 0.1e5)  # absolute stresses divided by pressure amplitude  #calculated originally from bottom to top - num idx guarantees that parameters from top to bottom are called
        sigma_num_xx = (sigma_xx[num_idx]/ 0.1e5)
        sigma_num_xy = (sigma_xy[num_idx]/ 0.1e5)

        f_abs_sigma_yy[num_idx] = sigma_num_yy - sigma_ana_yy     # saving with num_idx guarantees a saving from bottom to top (later it is plotted from bottom to top)
        f_abs_sigma_xx[num_idx] = sigma_num_xx - sigma_ana_xx
        f_abs_sigma_xy[num_idx] = sigma_num_xy - sigma_ana_xy

    return f_abs_sigma_yy,f_abs_sigma_xx, f_abs_sigma_xy

# %%
x = 0
y = np.linspace(0, 100, 1000)
y_rel = y / 100

fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(15, 15))

## Plotting analytical solution
for t in [2, 4, 6, 8, 10]:
    if t !=10:
        ax[0][0].plot(
            compute_pressure_and_stresses(t, x, y)[0],
            -y_rel,
            color=colors[t],
            alpha = 0.6,
        )
    else:
        ANALYTICAL, = ax[0][0].plot(
            compute_pressure_and_stresses(t, x, y)[0],
            -y_rel,
            color=colors[t],
            alpha = 0.6,
            label=f"analytical",
        )

ax[1][0].plot(
    compute_pressure_and_stresses(2, x, y)[1],
    -y_rel,
    color=colors[6],
    label="analytical, $\\sigma'_{xx}$ / $\\alpha\\tilde{p}$",
)
ax[1][0].plot(
    compute_pressure_and_stresses(2, x, y)[2],
    -y_rel,
    color=colors[2],
    label="analytical, $\\sigma'_{yy}$ / $\\alpha\\tilde{p}$",
)
ax[1][0].plot(
    compute_pressure_and_stresses(2, x, y)[3],
    -y_rel,
    color=colors[4],
    label="analytical, $\\sigma'_{xy}$ / $\\alpha\\tilde{p}$",
)

i_steps=[2,4,6,8,10]
meshes = [mesh_t2,mesh_t4,mesh_t6,mesh_t8,mesh_t10]
meshes_stress = [mesh_t2_stress]

for step, mesh in zip(i_steps, meshes):
    mask, line_points = nodes_along_line(mesh)
    pressure = get_pressure_sorted(mesh,mask,line_points)
    depth = get_depth_sorted_nodes(line_points)
    f_abs_pressure = compute_abs_and_rel_pressure_error(pressure, depth, step, x)

    if step != 10:
        ax[0][0].plot(
            pressure/0.1e5,
            (depth/100)-1,
            'o',
            color=colors[step]
            )
    else:
        ax[0][0].plot(
            pressure/0.1e5,
            (depth/100)-1,
            'o',
            color=colors[step],
            label=f"numerical - FEniCSx"
            )

    ax[0][0].set_xlabel("$p$ / $\\tilde{p}$")
    ax[0][0].set_ylabel("$y$ / L")
    ax[0][0].set_xlim(-1.1, 1.1)

    ax[0][1].plot(
        f_abs_pressure, 
        (depth / 100)-1, 
        'o',
        color=colors[step], 
        label=f"t = {step:.1f} s"
        )
    ax[0][1].set_xlabel("$\\Delta p$ / $\\tilde{p}$")
    ax[0][0].set_ylabel("$y$ / L")

# Add the first legend
first_patch = mpatches.Patch(color=colors[2], label='$t = 2.0$ s')
second_patch = mpatches.Patch(color=colors[4], label='$t = 4.0$ s')
third_patch = mpatches.Patch(color=colors[6], label='$t = 6.0$ s')
forth_patch = mpatches.Patch(color=colors[8], label='$t = 8.0$ s')
fifth_patch = mpatches.Patch(color=colors[10], label='$t = 10.0$ s')

first_legend = ax[0][0].legend(handles=[first_patch, second_patch, third_patch,forth_patch,fifth_patch], loc='lower left')

# Manually add the first legend back to the plot
ax[0][0].add_artist(first_legend)
ax[0][0].legend(bbox_to_anchor=(.6, 0.15))

for step, mesh_stress in zip(i_steps,meshes_stress):
    
    mask,line_points = nodes_along_line(mesh_stress)
    sigma_xx, sigma_yy, sigma_xy = get_stresses_sorted(mesh_stress,mask,line_points)
    depth =  get_depth_sorted_nodes(line_points)

    f_abs_sigma_yy,f_abs_sigma_xx,f_abs_sigma_xy= compute_abs_and_rel_stress_error(sigma_yy,sigma_xx,sigma_xy, depth, step, x)

    if step == 2:
        ax[1][0].plot(
            sigma_xx / 0.1e5,
            (depth / 100)-1,
            "o",
            color=colors[6],
            label="numerical, $\\sigma'_{xx}$ / $\\alpha\\tilde{p}$",
        )
        ax[1][0].plot(
            sigma_yy / 0.1e5,
            (depth / 100)-1,
            "o",
            color=colors[2],
            label="numercal, $\\sigma'_{yy}$ / $\\alpha\\tilde{p}$",
        )
        ax[1][0].plot(
            sigma_xy/0.1e5,
            (depth / 100)-1,
            "o",
            color=colors[4],
            label="numerical, $\\sigma'_{xy}$ / $\\alpha\\tilde{p}$",
        )
        ax[1][0].set_xlabel("$\\sigma$' / $\\alpha\\tilde{p}$")
        ax[1][0].set_ylabel("$y$ / L")
        ax[1][0].set_xlim(-1.1, 1.1)
        ax[1][0].text(.85, 0.05,"$t = 2.0$ s",
            transform=ax[1][0].transAxes,
            bbox=dict(
                facecolor='white', 
                edgecolor='black', 
                boxstyle='round,pad=0.5', 
                alpha=0.9))
        ax[1][0].legend()

        ax[1][1].plot(
            f_abs_sigma_yy,
            (depth / 100)-1,
            'o',
            color=colors[2],
            label="$\\Delta\\sigma'_{yy}$ / $\\alpha\\tilde{p}$",
        )
        ax[1][1].plot(
            f_abs_sigma_xx,
            (depth / 100)-1,
            'o',
            color=colors[6],
            label="$\\Delta\\sigma'_{xx}$ / $\\alpha\\tilde{p}$",
        )
        ax[1][1].plot(
            f_abs_sigma_xy,
            (depth / 100)-1,
            'o',
            color=colors[4],
            label="$\\Delta\\sigma'_{xy}$ / $\\alpha\\tilde{p}$",
        )
        ax[1][1].set_xlabel("$\\Delta\\sigma'$ / $\\alpha\\tilde{p}$")

plt.savefig('./'+name_result_folder+'/comparison_analytical_fenicsx_p_sigma_error.pdf')


