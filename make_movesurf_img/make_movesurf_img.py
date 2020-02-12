# flake8: noqa

import geonomics as gnx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

# define number of individuals to plot tracks for, and number of timesteps for
# tracks
n_individs = 20
n_timesteps = 150

# make figure
fig = plt.figure(figsize=(9.25, 4.5))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.065])

# make model
mod = gnx.make_model(('/home/drew/Desktop/stuff/berk/research/projects/'
                      'sim/methods_paper/movesurf_img/movesurf_img_params.py'))
# plot the movement_surface
ax1 = plt.subplot(gs[0])
mod.plot_movement_surface(0, 'chist', ticks=False)
ax1.set_title('von Mises mixture histograms', fontsize=30)

# plot tracks
ax2 = plt.subplot(gs[1])
im = plt.pcolormesh(np.linspace(0, 7, 8), np.linspace(0, 7, 8),
                    mod.land[0].rast, cmap='plasma')
gnx.help.plot_movement(mod.comm[0], mod.land, n_timesteps,
                       0, mod.params, subset_spp=n_individs-1,
                       ticks=False, color='gray', color_by_individ=False,
                       increasing_linewidth=False, alpha=0.5,
                       include_start_points=False)
gnx.help.plot_movement(mod.comm[0], mod.land, n_timesteps,
                       0, mod.params, subset_spp=1, ticks=False,
                       increasing_linewidth=False, alpha=0.7, color='black',
                       include_start_points=False)
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im, cax=cax)
cbar.set_label('conductance', rotation=270, labelpad=25, y=0.5, fontsize=24)
#ax2.set_title('Sample movement tracks\nfor %i individuals' % n_individs)
ax2.set_title('Sample movement tracks', fontsize=30)

fig.tight_layout()
plt.show()
fig.savefig(('/home/drew/Desktop/stuff/berk/research/projects/sim/'
             'methods_paper/img/final/move_surf.pdf'),
            format='pdf', dpi=1000)
