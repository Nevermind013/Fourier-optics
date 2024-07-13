import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

params = {'font.family': 'Times New Roman', 'legend.fontsize': 18, 'figure.figsize': (9.5, 6), 'axes.labelsize': 18,
          'axes.titlesize': 20, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'figure.titlesize': 20}
plt.rcParams.update(params)
plt.style.use(['science', 'notebook'])
sf = ScalarFormatter()


def screen(xv, yv, zv, title, sfx, sfy):
    """
    visualise the intensity profile zv in the oxy plane

    :param xv: ox grid
    :param yv: oy grid
    :param zv: 2d-array described as zv=func(xv, yv)
    :param title: title for the image
    :param sfx: scalarformatter for ox
    :param sfy: scalarformatter for oy
    :return: None
    """

    fig, ax = plt.subplots(figsize=(6, 5), edgecolor='black', linewidth=3, frameon=True)
    im = ax.pcolormesh(xv, yv, zv, cmap='inferno')
    ax.set_title(title)
    ax.set_xlabel('$x$ [mm]')
    ax.set_ylabel('$y$ [mm]')

    ax.xaxis.set_major_formatter(sfx)
    ax.yaxis.set_major_formatter(sfy)
    fig.colorbar(im)
    return None


def relative_error(intensity, target):
    """
    count relative error between target intensity T

    :param intensity: calculated intensity profile
    :param target: target intensity profile
    :return: relative error
    """

    region = (target > np.amax(target) / 2).astype(float)
    region_num = int(np.sum(region))

    norm_target = np.sum(target * region)
    norm_intensity = np.sum(intensity * region)

    target_normalised = target * region / norm_target
    intensity_normalised = intensity * region / norm_intensity

    error = np.sqrt((1 / region_num) * np.sum(np.power(target_normalised*region - intensity_normalised*region, 2) /
                    np.power(target_normalised*region+1e-22, 2)))
    return error



