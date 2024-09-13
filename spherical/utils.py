import os, math, numpy as np, matplotlib.pyplot as plt
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


def air_mass(zenith_angle: float) -> float:
    """
    calculate the air mass for a given zenith angle
    """
    if zenith_angle > 89.0:
        return 40.0 # approximate max air mass near the horizon
    
    return 1.0 / (math.cos(math.radians(zenith_angle)) + 0.15 * (93.885 - zenith_angle) ** -1.253) # kasten-young model


def refraction_correction(zenith_angle: float) -> float:
    """
    apply a basic atmospheric refraction correction for light near the horizon
    """
    if zenith_angle < 90.0:
        z = math.radians(zenith_angle)
        # approximate refraction formula for standard atmosphere
        return (1.02 / math.tan(z + (10.3 / (z + 5.11)))) / 60.0 # returns refraction in degrees
    
    return 0.0


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


def gaussian_kernel(size: int, sigma: float) -> np.ndarray:
    """
    create a 1D gaussian kernel using the specified size and standard deviation (sigma)

    args:
        size: size of the gaussian kernel (typically an odd number)
        sigma: standard deviation of the gaussian distribution

    returns:
        1D numpy array representing the gaussian kernel
    """
    # create an array of values [-size // 2, size // 2]
    x = np.arange(-size // 2 + 1, size // 2 + 1)

    # compute the 1D gaussian distribution for these values
    kernel = np.exp(-x ** 2 / (2 * sigma ** 2))

    # normalize the kernel to ensure the sum is 1
    return kernel / np.sum(kernel)


def apply_gaussian_smoothing(data: np.ndarray, sigma: float) -> np.ndarray:
    """
    apply a gaussian smoothing filter to the input data

    args:
        data: 1D numpy array of intensity values to smooth
        sigma: standard deviation of the gaussian kernel

    returns:
        smoothed 1D numpy array
    """
    # choose the size of the gaussian kernel (3 * sigma is a good rule of thumb)
    kernel_size = int(3 * sigma) * 2 + 1

    # generate the gaussian kernel
    kernel = gaussian_kernel(kernel_size, sigma)

    # convolve the data with the kernel using 'same' mode to keep the same size
    smoothed_data = np.convolve(data, kernel, mode='same')

    return smoothed_data
