from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import numpy as np
import math
import argparse
from surface import Surface, create_surfaces
from quantum import get_probability_fn, get_density_fn, get_radius, a0
from point import Point, cart2sph, sph2cart

# Z is the atomic number
# n, l, and m are quantum numbers
# r, t, and p are the spherical coordinates
#     with ranges r:[0,inf), t:[0,pi], p:[0,2*pi] 
# x, y, z are the Cartesian coordinates

def get_spherical_values(fn, max_r, res):
    # Returns matrix of fn vals and axis vals of spherical coordinates
    vals = np.empty((res, 2*res, res))
    for i, t in enumerate(np.linspace(0, math.pi, vals.shape[0])):
        for j, p in enumerate(np.linspace(0, 2*math.pi, vals.shape[1])):
            for k, r in enumerate(np.linspace(0, max_r, vals.shape[2])):
                vals[i,j,k] = fn(r, t, p)
    return vals

def get_threshold(data, percentage):
    # Return the value in data for the percentage
    s_data = np.sort(data.flatten())
    limit = np.sum(s_data) * percentage
    total = 0
    for threshold in s_data[::-1]:
        total += threshold
        if total >= limit:
            break
    return threshold

def graph_2d(n, l, m, ax, intersect=lambda x,y: (x,y,0), res=128):
    # intersect is the function for the intersection surface
    max_r = get_radius(get_probability_fn(n, l, m))
    # Populate data with probability densities
    fn = get_density_fn(n, l, m)
    data = np.empty((res,res))
    for i, x in enumerate(np.linspace(-1*max_r, max_r, res)):
        for j, y in enumerate(np.linspace(-1*max_r, max_r, res)):
            data[j, i] = fn(*cart2sph(*intersect(x, y)))
    # Scale probability densities to probabilities per pixel
    data /= np.sum(data)
    # Draw the graph
    ax.imshow(data, cmap="hot", extent=[-1*max_r/a0, max_r/a0, -1*max_r/a0, max_r/a0])

def graph_3d(n, l, m, ax, percentages=[0.25, 0.5, 0.75], res=32):
    # Plot the surfaces
    P = get_probability_fn(n, l, m)
    D = get_density_fn(n, l, m)
    max_r = get_radius(P)
    vals = get_spherical_values(D, max_r, res)
    for i, percentage in enumerate(sorted(percentages)):
        threshold = get_threshold(vals, percentage)
        fn = lambda r, t, p: D(r, t, p) > threshold
        for surface in create_surfaces(fn, max_r, res):
            ax.plot_wireframe(surface.x(), surface.y(), surface.z(), cstride=res/32, rstride=res/32, alpha=1-0.6/(len(percentages)-1)*i, color=cm.hot(1-percentage))
    # Configure the plot
    # Black background
    ax.set_axis_bgcolor((0,0,0))
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')
    # White axis labels
    ax.xaxis.label.set_color('w')
    ax.yaxis.label.set_color('w')
    ax.zaxis.label.set_color('w')
    # White axis
    ax.tick_params(axis='x', colors='w')
    ax.tick_params(axis='y', colors='w')
    ax.tick_params(axis='z', colors='w')
    # Black axis background
    ax.w_xaxis.set_pane_color((0,0,0,1))
    ax.w_yaxis.set_pane_color((0,0,0,1))
    ax.w_zaxis.set_pane_color((0,0,0,1))
    # Axis units are a0
    ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/a0))
    ax.xaxis.set_major_formatter(ticks)
    ax.yaxis.set_major_formatter(ticks)
    ax.zaxis.set_major_formatter(ticks)
    # Make axis equal
    lims = zip(ax.get_xlim(), ax.get_ylim(), ax.get_zlim())
    min_lim, max_lim = min(lims[0]), max(lims[1])
    ax.set_xlim([min_lim, max_lim])
    ax.set_ylim([min_lim, max_lim])
    ax.set_zlim([min_lim, max_lim])
    # Make ticks in whole numbers of a0
    ticks = np.linspace(math.ceil(min_lim/a0), math.floor(max_lim/a0), 5)*a0
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_zticks(ticks)

def show_graphs(n, l, m):
    print 'Calculating graphs, this can take a while...'
    fig = plt.figure(num=None, figsize=(8, 8), dpi=80, facecolor='w', edgecolor='k')
    fig.canvas.set_window_title('Hydrogen probability density graphs for n={}, l={}, m={} (Units in Bohr radius)'.format(n, l, m))
    ax1 = plt.subplot2grid((4, 3), (0, 0))
    ax2 = plt.subplot2grid((4, 3), (0, 1))
    ax3 = plt.subplot2grid((4, 3), (0, 2))
    ax4 = plt.subplot2grid((4, 3), (1, 0), colspan=3, rowspan=3, projection='3d')
    graph_2d(n, l, m, ax1, lambda x,y: (x,y,0))
    ax1.set_xlabel('X-axis')
    ax1.set_ylabel('Y-axis')
    graph_2d(n, l, m, ax2, lambda x,z: (x,0,z))
    ax2.set_xlabel('X-axis')
    ax2.set_ylabel('Z-axis')
    graph_2d(n, l, m, ax3, lambda y,z: (0,y,z))
    ax3.set_xlabel('Y-axis')
    ax3.set_ylabel('Z-axis')
    graph_3d(n, l, m, ax4)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create and display hydrogen probability density graphs.')
    parser.add_argument('n', metavar='n', type=int, choices=[1,2,3,4,5],
                    help='Principal quantum number')
    parser.add_argument('l', metavar='l', type=int, choices=[0,1,2,3,4],
                        help='Azimuthal quantum number')
    parser.add_argument('m', metavar='m', type=int, choices=[-4,-3,-2,-1,0,1,2,3,4],
                        help='Magnetic quantum number')
    args = parser.parse_args()
    if args.l >= args.n:
        parser.error('The azimuthal qunatum number, l, must be less than the principal quantum number, n.')
    if abs(args.m) > args.l:
        parser.error('The absolute value of the magnetic quantum number, m, cannot be greater than the azimuthal quantum number, l.')
    show_graphs(args.n,args.l,args.m)
 