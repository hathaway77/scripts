from pythtb import *
import matplotlib.pyplot as plt
from matplotlib import cm
# lattice vectors and orbital positions
lat = [[1.0, 0.0], [-0.5, np.sqrt(3.0)/2.0]]
orb = [[2./3., 0.], [0., 0.], [0., 1./3.]]
# make two dimensional tight-binding model
my_model = tb_model(2, 2, lat, orb)
# set model parameters
delta = 0.2
tnn = -1.0
tnnn = 1/3*tnn
# set on-site energies
my_model.set_onsite([0., 0., 0.])
# set hoppings (one for each connected pair of orbitals)
# (amplitude, i, j, [lattice vector to cell containing j])
my_model.set_hop(tnn, 0, 1, [1, 0])
my_model.set_hop(tnn, 0, 2, [1, 0])
my_model.set_hop(tnn, 1, 2, [0, 0])

my_model.set_hop(tnnn, 0, 1, [0, 0])
my_model.set_hop(tnnn, 0, 2, [0, -1])
my_model.set_hop(tnnn, 1, 2, [0, -1])
# print tight-binding model
my_model.display()
# generate k-point path and labels
path=[[0.0,0.0],[0.0,0.5],[0.33333,0.33333],[0.0,0.0]]
label=(r'$\Gamma $',r'$M$', r'$K$', r'$\Gamma $')
(k_vec,k_dist,k_node)=my_model.k_path(path,301,report=False)
#print(k_dist)
#print(k_node)
evals=my_model.solve_all(k_vec)

# First make a figure object
fig, ax = plt.subplots()
# specify horizontal axis details
ax.set_xlim(k_node[0],k_node[-1])
ax.set_xticks(k_node)
ax.set_xticklabels(label)
for n in range(len(k_node)):
    ax.axvline(x=k_node[n], linewidth=0.5, color='k')

# plot bands
for n in range(3):
    ax.plot(k_dist,evals[n])
# put title
ax.set_title("Breathing_Kagome band structure")
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy\nE(tnn)")
# make an figure of a plot
fig.tight_layout()
plt.savefig('Kagome_2D.png')

#3D
num = 101
X = np.linspace(-5, 5, num)
Y = np.linspace(-5, 5, num)

evals_a_3d = []
evals_b_3d = []
evals_c_3d = []

for i in X:
    for j in Y:
        k_3d = [i*2.0/np.sqrt(3.0)/(2*np.pi), (j-i/np.sqrt(3.0))/(2*np.pi)]
        eval=my_model.solve_one(k_3d)
        evals_a_3d.append(eval[0])
        evals_b_3d.append(eval[1])
        evals_c_3d.append(eval[2])

evals_a_3d = np.array(evals_a_3d)
evals_a_3d.resize(num, num)
evals_b_3d = np.array(evals_b_3d)
evals_b_3d.resize(num, num)
evals_c_3d = np.array(evals_c_3d)
evals_c_3d.resize(num, num)
#网格化
X, Y = np.meshgrid(X,Y)
a = [np.pi/np.sqrt(3.0), 2*np.pi/np.sqrt(3.0), np.pi/np.sqrt(3.0), -np.pi/np.sqrt(3.0), -2*np.pi/np.sqrt(3.0), -np.pi/np.sqrt(3.0), np.pi/np.sqrt(3.0)]
b = [np.pi, 0, -np.pi, -np.pi, 0, np.pi, np.pi]
c = [-1, -1, -1, -1, -1, -1, -1]
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
#plt.gca().set_box_aspect((1.65, 0.952, 1.5))
#ax.plot(b, a, c, linewidth=2, color='black')
surf = ax.plot_surface(X, Y, evals_a_3d, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
surf2 = ax.plot_surface(X, Y, evals_b_3d, rstride=1, cstride=1, cmap=cm.RdYlGn,
        linewidth=0, antialiased=False)
surf3 = ax.plot_surface(X, Y, evals_c_3d, rstride=1, cstride=1, cmap=cm.cool,
        linewidth=0, antialiased=False)

#ax.contour(X, Y, evals_a_3d, zdir = 'z', offset = -4, cmap = plt.get_cmap('rainbow'))
ax.set_xlim3d(-6, 6)
ax.set_ylim3d(-6, 6)
ax.set_zlabel("Band energy\nE(tnn)")
#ax.set_zlim3d(-4, 2)
#fig.colorbar(surf, shrink=0.2, aspect=10)
ax.set_title("Breathing_Kagome 3D band structure")
ax.view_init(elev=15,azim=120)
#####################
plt.savefig("Kagome_3D.png",dpi=600)

