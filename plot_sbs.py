from astropy.io import fits as pf
import numpy as np
import matplotlib.pyplot as plt
import math as m
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib.ticker as tkr


def load_data(files, scale):
    """ 
    Load FITS files and convert the units with a scale factor
    """
    data = {}
    # Load files
    for file_name in files:
        # Load the FITS
        raw_data = pf.getdata(files[file_name])

        # Minimal processing and scaling
        raw_data = np.nan_to_num(raw_data[0,0,:,:]) * scale
        data[file_name] = raw_data

    return data

data_paths = {
    'lp349': '349_concat_clean.fits', 
    'lsr': 'cal_final_tclean.fits',
    'nltt': 'NLTT_clean.fits',
}

def add_to_plot(
        plot_data,
        plot_index,
        resolution,
        name,
        beam_x,
        beam_y,
        beam_pa,
):
    """
    Plot a single panel of the continuum plots
    """

    # Create an ellipse representing the synthetic beam
    ellipse = Ellipse((-8,-8), beam_y, beam_x, angle = beam_pa, color = 'black')
    ax = fig.add_subplot(1,3,plot_index, aspect = "equal")
    ax.add_patch(ellipse)
    minorLocator = MultipleLocator(5)

    # Create the contour data
    levels = np.linspace(np.min(plot_data), np.max(plot_data), num=200)
    lengthx = resolution * plot_data.shape[0] 
    lengthy = resolution * plot_data.shape[1] 
    dat =ax.contourf(
        plot_data, 
        origin = 'lower', 
        levels = levels, 
        extent = [-lengthx/2., lengthx/2., -lengthy/2., lengthy/2.], 
        cmap='coolwarm'
    )

    # Additional plot formatting
    plt.locator_params(axis='y', nbins=4)
    plt.locator_params(axis='x', nbins=4)
    plt.xlabel(r'$\Delta \alpha$ ["]', fontsize = 10)
    plt.ylabel(r'$\Delta \delta$ ["]', fontsize = 10)
    plt.xlim(-10,10)
    plt.ylim(-10,10)
    plt.title(name, fontsize = 15)
    ax.set_xticklabels([10, 5, 0, -5, -10]) # Reverse x-axis convention
    ax.xaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_label_coords(-0.07,0.5)
    ax.tick_params(axis='both', which='major', labelsize=10)

    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(dat, ax=ax, cax=cax, format=tkr.FormatStrFormatter('%.0f'))
    cbar.set_label(r'$\rm \mu Jy  \, beam^{-1}$', labelpad = -4, fontsize = 8)


# Load the data
data = load_data(files = data_paths, scale = 1000000)

# Define some plot parameters
plot_data = {
    'lp349': {
        'name': 'LP 349-35',
        'beam_x': 0.823,
        'beam_y': 0.621,
        'beam_pa': 28.3,
        'resolution': 0.099
    },
    'lsr': {
        'name': 'LSR J1835+3259',
        'beam_x': 2.853,
        'beam_y': 1.185,
        'beam_pa': -41.8,
        'resolution': 0.11
    },
    'nltt': {
        'name': 'NLTT 33370',
        'beam_x': 1.567,
        'beam_y': 1.037,
        'beam_pa': -44.4,
        'resolution': 0.099
    }
}


# Loop through the data to create plots
counter = 1
fig = plt.figure(figsize=(12,4))
for data_id in data:
    add_to_plot(
        plot_data=data[data_id], 
        plot_index=counter, 
        resolution=plot_data[data_id]['resolution'], 
        name=plot_data[data_id]['name'], 
        beam_x=plot_data[data_id]['beam_x'],
        beam_y=plot_data[data_id]['beam_y'],
        beam_pa=plot_data[data_id]['beam_pa']
        )
    counter+=1


# Save
fig.tight_layout(pad=0.1)
plt.savefig("continuum_plots.png")
plt.show()
