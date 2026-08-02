import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema as extrema

CM = 1 / 2.54  # cm -> inches conversion factor

def make_stacked_axes(n_panels, panel_height_cm=5.0, gap_cm=0.3,
                       fig_width_cm=16, left=0.15, right=0.95,
                       top_margin_cm=0.3, bottom_margin_cm=1.5):
    """
    Create n_panels stacked subplots (shared x-axis), each with the
    SAME fixed height in cm, regardless of how many panels there are.
    The figure height adapts automatically.

    Parameters
    ----------
    n_panels : int
        Number of stacked subplots.
    panel_height_cm : float
        Height of each individual panel, in cm.
    gap_cm : float
        Vertical gap between panels, in cm.
    fig_width_cm : float
        Total figure width, in cm.
    left, right : float
        Left/right margins as fractions of figure width (0 to 1).
    top_margin_cm, bottom_margin_cm : float
        Extra space (in cm) reserved at the very top (e.g. for a title)
        and very bottom (e.g. for x-axis tick labels).

    Returns
    -------
    fig : Figure
    axes : list of Axes, ordered top to bottom
    """
    # total figure height = margins + all panels + all gaps between them
    total_height_cm = (
        top_margin_cm
        + n_panels * panel_height_cm
        + (n_panels - 1) * gap_cm
        + bottom_margin_cm
    )

    fig = plt.figure(figsize=(fig_width_cm * CM, total_height_cm * CM))

    frac_panel = panel_height_cm / total_height_cm
    frac_gap = gap_cm / total_height_cm
    frac_top_margin = top_margin_cm / total_height_cm

    axes = []
    top = 1.0 - frac_top_margin  # start below the top margin

    for i in range(n_panels):
        bottom = top - frac_panel
        ax = fig.add_axes((left, bottom, right - left, frac_panel))
        axes.append(ax)
        top = bottom - frac_gap

    # share x-axis and hide tick labels on all but the bottom panel
    for ax in axes[:-1]:
        ax.sharex(axes[-1])
        ax.tick_params(labelbottom=False)

    return fig, axes



def gapWidth(x_data,y_data):
    n_points = len(y_data)

    left_x = x_data[:n_points//2]
    right_x = x_data[n_points//2:]

    left_y = y_data[:n_points//2]
    right_y = y_data[n_points//2:]

    left_index = extrema(np.array(left_y),np.greater)
    right_index = extrema(np.array(right_y),np.greater)

    max_left = left_index[0][-1]
    max_right = right_index[0][0]

    max_x_left = left_x[max_left]
    max_x_right = right_x[max_right]

    return max_x_left, max_x_right
