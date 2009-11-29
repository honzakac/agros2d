# model
newdocument("Rotating induction heating", "planar", "magnetostatic", 1, 2, "disabled", 3, 0, 0, "steadystate", 1, 1, 0)

# coil
Jext = 3e7

# heated part
sigma = 33e6
n = 1500
p = 1
f = p*n/60
omega = 2*pi*f

# boundaries
addboundary("A = 0", "magnetostatic_vector_potential", 00)

# materials
addmaterial("Al", 0, 1, 0, 0, sigma, 0, 0, omega)
addmaterial("Air", 0, 1, 0, 0, 0, 0, 0, 0)
addmaterial("Cu+", Jext, 1, 0, 0, 0, 0, 0, 0)
addmaterial("Cu-", -Jext, 1, 0, 0, 0, 0, 0, 0)

# edges
addedge(0.03, -0.0519615, 0.03, 0.0519615, 120, "none")
addedge(0.03, 0.0519615, 0.04, 0.069282, 0, "none")
addedge(0.04, -0.069282, 0.04, 0.069282, 120, "none")
addedge(0.04, -0.069282, 0.03, -0.0519615, 0, "none")
addedge(-0.03, -0.0519615, -0.04, -0.069282, 0, "none")
addedge(-0.03, 0.0519615, -0.04, 0.069282, 0, "none")
addedge(-0.04, 0.069282, -0.04, -0.069282, 120, "none")
addedge(-0.03, 0.0519615, -0.03, -0.0519615, 120, "none")
addedge(0.05, 0, 0, 0.05, 90, "none")
addedge(0, 0.05, -0.05, 0, 90, "none")
addedge(-0.05, 0, 0, -0.05, 90, "none")
addedge(0, -0.05, 0.05, 0, 90, "none")
addedge(0.25, 0, 0, 0.25, 90, "A = 0")
addedge(0, 0.25, -0.25, 0, 90, "A = 0")
addedge(-0.25, 0, 0, -0.25, 90, "A = 0")
addedge(0, -0.25, 0.25, 0, 90, "A = 0")

# labels
addlabel(0.0899971, 0.155764, 0, "Air")
addlabel(-0.0184609, 0.0207686, 0, "Al")
addlabel(0.0576904, 0.0403833, 0, "Cu+")
addlabel(-0.0588442, 0.0380757, 0, "Cu-")

zoombestfit()
solve()

N = 30
step = 0.05/N

def graph(angle):
	r = []
	Jvel = []
	for i in range(1, N+1):
		R = (i-1)*step
		point = pointresult(R*cos(angle/180.0*pi), R*sin(angle/180.0*pi))
		r.append(R)
		Jvel.append(point["Jvel"])
	return r, Jvel


r00, Jvel00 = graph(00)
r30, Jvel30 = graph(30)
r60, Jvel60 = graph(60)
r90, Jvel90 = graph(90)

# plot chart
import pylab

pylab.clf()
# pylab.plot(r00, Jvel00, "k-", r30, Jvel30, "k+-", r60, Jvel60, "k--", r90, Jvel90, "ko-")
pylab.plot(r00, Jvel00, r30, Jvel30, r60, Jvel60, r90, Jvel90)
pylab.grid(1)
pylab.xlabel("r (m)")
pylab.ylabel("J (A/m2)")
pylab.legend(("00 deg.", "30 deg.", "60 deg.", "90 deg."), "upper left")
pylab.show()