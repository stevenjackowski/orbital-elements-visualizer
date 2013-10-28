__author__ = 'sajackow'

from Tkinter import Tk, DoubleVar, IntVar
from ttk import Frame, Label, Scale
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from mpl_toolkits.mplot3d import Axes3D, proj3d
from matplotlib.patches import FancyArrowPatch
from orbit import *


# Used for drawing arrows
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


class OrbitUI(Frame):
    def __init__(self, parent):
        Frame.__init__(self, parent)

        self.parent = parent
        self.initUI()
        self.pack()

    def initUI(self):
        self.parent.title("Orbit Visualizer")

        # tk variables corresponding to CoE's
        self.ecc_var = IntVar()
        self.ecc_var.set(0)
        self.inc_var = IntVar()
        self.inc_var.set(0)
        self.raan_var = IntVar()
        self.raan_var.set(0)
        self.aop_var = IntVar()
        self.aop_var.set(0)
        self.a_var = IntVar()
        self.a_var.set(7000)
        self.ta_var = IntVar()
        self.ta_var.set(0)
        self.ecc_actual = DoubleVar()
        # use ecc_actual to make the resolution of ecc = 0.01
        self.ecc_actual.set(self.ecc_var.get()/100.0)
        self.initOrbit()

        # eccentricity
        Label(self, text="e").grid(row=1)
        ecc_scale = Scale(self, from_=0, to=80, variable=self.ecc_var,
                          command=self.updateOrbit)
        ecc_scale.grid(row=1, column=1, columnspan=2)


        # inclination
        Label(self, text="i").grid(row=2)
        inc_scale = Scale(self, from_=-180, to=180, variable=self.inc_var,
                          command=self.updateOrbit)
        inc_scale.grid(row=2, column=1, columnspan=2)

        # RAAN
        Label(self, text=u"\u03A9").grid(row=3)
        raan_scale = Scale(self, from_=0, to=360, variable=self.raan_var,
                           command=self.updateOrbit)
        raan_scale.grid(row=3, column=1, columnspan=2)

        # argument of perigee
        Label(self, text=u"\u03c9").grid(row=4)
        aop_scale = Scale(self, from_=0, to=360, variable=self.aop_var,
                          command=self.updateOrbit)
        aop_scale.grid(row=4, column=1, columnspan=2)

        # semi-major axis
        Label(self, text="a").grid(row=5)
        a_scale = Scale(self, from_=6378, to=50000, variable=self.a_var,
                        command=self.updateOrbit)
        a_scale.grid(row=5, column=1, columnspan=2)

        # true anomaly
        Label(self, text=u"M\u02f3").grid(row=6)
        ta_scale = Scale(self, from_=0, to=6, variable=self.ta_var,
                         command=self.updateOrbit)
        ta_scale.grid(row=6, column=1, columnspan=2)

    def initOrbit(self):
        self.orbit = Orbit(self.ecc_actual.get(), self.inc_var.get(),
                           self.raan_var.get(), self.aop_var.get(),
                           self.a_var.get())

        # set up axis and labels
        self.fig = plt.figure(figsize=plt.figaspect(0.9))
        self.fig.patch.set_facecolor('white')
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_xlim3d([-self.a_var.get(), self.a_var.get()])
        self.ax.set_ylim3d([-self.a_var.get(), self.a_var.get()])
        self.ax.set_zlim3d([-self.a_var.get(), self.a_var.get()])
        self.ax.set_xlabel("x (km)")
        self.ax.set_ylabel("y (km)")
        self.ax.set_zlabel("z (km)")

        r_vecs = self.orbit.get_rvectors()
        self.orbit_data, = self.ax.plot(r_vecs[0, :], r_vecs[1, :], r_vecs[2,
                                                                   :],'r')
        # plot the Earth as a sphere
        u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:20j]
        x = np.cos(u)*np.sin(v)*6378
        y = np.sin(u)*np.sin(v)*6378
        z = np.cos(v)*6378
        self.ax.plot_wireframe(x, y, z, color="b")

        # #
        # add arrows and annotations
        # #
        # reference direction
        ref_arrow = Arrow3D([0, 6500], [0, 0], [0, 0], mutation_scale=20,
                            lw=1, arrowstyle="-|>", color="k")
        self.ax.add_artist(ref_arrow)

        # create tk canvas object
        self.orbit_plot = FigureCanvasTkAgg(self.fig, master=self)
        self.orbit_plot.show()
        self.orbit_plot.get_tk_widget().grid(row=0, column=5, rowspan=8,
                                             columnspan=8)

    def updateOrbit(self, x=0):
        self.ecc_actual.set(self.ecc_var.get()/100.0)
        self.orbit.update(self.ecc_actual.get(), self.inc_var.get(),
                           self.raan_var.get(), self.aop_var.get(),
                           self.a_var.get())
        r_vecs = self.orbit.get_rvectors()

        # update the plot data
        self.orbit_data.set_data(r_vecs[0, :], r_vecs[1, :])
        self.orbit_data.set_3d_properties(r_vecs[2, :])
        self.ax.set_xlim3d([-self.a_var.get(), self.a_var.get()])
        self.ax.set_ylim3d([-self.a_var.get(), self.a_var.get()])
        self.ax.set_zlim3d([-self.a_var.get(), self.a_var.get()])
        self.ax.auto_scale_xyz(r_vecs[0, :], r_vecs[1, :], r_vecs[2, :])
        self.orbit_plot.draw()

def main():
    root = Tk()
    app = OrbitUI(root)
    root.mainloop()

if __name__ == '__main__':
    main()