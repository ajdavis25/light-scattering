import numpy as np


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
