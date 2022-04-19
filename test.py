import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def main():
    # Generate some random data
    nx, ny = 100, 100
    data = np.random.random((ny,nx))

    # Define a circle in the center of the data with a radius of 20 pixels
    radius = 20
    center_x = nx // 2
    center_y = ny // 2

    plot_masked(data, center_x, center_y, radius)
    plot_clipped(data, center_x, center_y, radius)
    plt.show()

def plot_masked(data, center_x, center_y, radius):
    """Plots the image masked outside of a circle using masked arrays"""
    # Calculate the distance from the center of the circle
    ny, nx = data.shape
    ix, iy = np.meshgrid(np.arange(nx), np.arange(ny))
    distance = np.sqrt((ix - center_x)**2 + (iy - center_y)**2)

    # Mask portions of the data array outside of the circle
    data = np.ma.masked_where(distance > radius, data)

    # Plot
    plt.figure()
    plt.imshow(data)
    plt.title('Masked Array')

def plot_clipped(data, center_x, center_y, radius):
    """Plots the image clipped outside of a circle by using a clip path"""
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Make a circle
    circ = patches.Circle((center_x, center_y), radius, facecolor='none')
    ax.add_patch(circ) # Plot the outline

    # Plot the clipped image
    im = ax.imshow(data, clip_path=circ, clip_on=True)

    plt.title('Clipped Array')

main()

