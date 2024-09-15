import os, numpy as np, matplotlib.pyplot as plt
from typing import Dict, List, Tuple
from matplotlib.colors import LinearSegmentedColormap


def save_plot(fig, save_path: str, save_name: str) -> None:
    """
    saves a plot figure to a specified location and optionally displays it

    args:
        fig (matplotlib.figure.Figure): the Matplotlib figure object to save
        save_path (str): the path to the directory where the plot will be saved
        save_name (str): the filename for the saved plot
    """
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    plt.savefig(os.path.join(save_path, save_name))
    plt.show()


def create_twilight_colormap() -> LinearSegmentedColormap:
    """
    create a custom colormap to represent the twilight sky

    returns:
        a LinearSegmentedColormap object
    """
    # define key colors and positions in the gradient (0 = far from sun, 1 = near the sun)
    cdict: Dict[str, List[Tuple[float, float, float]]] = {
        'red':   [(0.0, 0.0, 0.0),   # black (night sky)
                  (0.2, 0.2, 0.2),   # dark blue
                  (0.5, 0.5, 0.5),   # purple
                  (0.7, 0.9, 0.9),   # red/orange (sunset colors)
                  (1.0, 1.0, 1.0)],  # yellow (bright sunlight near horizon)
        
        'green': [(0.0, 0.0, 0.0),   # black
                  (0.2, 0.1, 0.1),   # dark blue
                  (0.5, 0.0, 0.0),   # purple
                  (0.7, 0.5, 0.5),   # red/orange
                  (1.0, 1.0, 1.0)],  # yellow
        
        'blue':  [(0.0, 0.0, 0.0),   # black
                  (0.2, 0.5, 0.5),   # dark blue
                  (0.5, 0.5, 0.5),   # purple
                  (0.7, 0.0, 0.0),   # red/orange
                  (1.0, 0.0, 0.0)]   # yellow
    }

    # create the colormap object
    twilight_colormap = LinearSegmentedColormap('twilight_sky', cdict)

    return twilight_colormap
