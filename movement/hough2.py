# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

# Create synthetic noisy data (a line with noise)
np.random.seed(0)
n_points = 100

# Generate points along a line y = mx + c, where m = 2, c = 10
x = np.linspace(0, 200, n_points)
y = 2 * x + 10

# Add Gaussian noise
x_noise = x + np.random.normal(0, 20, n_points)
y_noise = y + np.random.normal(0, 20, n_points)

# Add extra noise 
x_noise = np.concatenate((x_noise, np.random.normal(0, max(x_noise)/2, 200)))
y_noise = np.concatenate((y_noise, np.random.normal(0, max(y_noise)/2, 200)))

# Plot the noisy data
plt.figure()
plt.scatter(x_noise, y_noise, color='b', label='Noisy data')
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.show()





# Hough Transform Implementation
def hough_transform(x_points, y_points, theta_res=1, rho_res=1):
    # Convert angles to radians
    theta_range = np.deg2rad(np.arange(-90, 90, theta_res))
    
    # Create an empty accumulator for rho and theta values
    max_rho = int(np.sqrt(x_points.max()**2 + y_points.max()**2))
    rhos = np.arange(-max_rho, max_rho, rho_res)
    accumulator = np.zeros((len(rhos), len(theta_range)))

    # Hough Transform: vote in the accumulator for each point
    for x, y in zip(x_points, y_points):
        for theta_index, theta in enumerate(theta_range):
            rho = x * np.cos(theta) + y * np.sin(theta)
            rho_index = np.argmin(np.abs(rhos - rho))  # Find the closest rho bin
            accumulator[rho_index, theta_index] += 1   # Increment the accumulator

    return accumulator, theta_range, rhos

# Run the Hough Transform
accumulator, theta_range, rhos = hough_transform(x_noise, y_noise, 1, 10)

# Display the Hough accumulator
plt.figure()
plt.imshow(accumulator, extent=[np.rad2deg(theta_range.min()), np.rad2deg(theta_range.max()), rhos.min(), rhos.max()], aspect='auto')
plt.xlabel('Theta (degrees)')
plt.ylabel('Rho')
plt.title('Hough Space')
plt.colorbar(label='Votes')
plt.show()




# Find the peaks in the accumulator
def find_hough_peaks(accumulator, num_peaks=1, threshold=100):
    peaks = []
    for _ in range(num_peaks):
        rho_idx, theta_idx = np.unravel_index(np.argmax(accumulator), accumulator.shape)
        if accumulator[rho_idx, theta_idx] >= threshold:
            peaks.append((rho_idx, theta_idx))
            accumulator[rho_idx, theta_idx] = 0  # Zero out this peak
    return peaks

def find_hough_peak(accumulator):
    x, y = np.where(accumulator == np.max(accumulator))
    return (x[0], y[0])
    


# Detect the strongest peaks
peaks = find_hough_peaks(accumulator, num_peaks=1, threshold=10)

print('Peaks: ', peaks)


# Visualize the detected peaks in Hough space
for rho_idx, theta_idx in peaks:
    plt.scatter(np.rad2deg(theta_range[theta_idx]), rhos[rho_idx], color='red')

plt.figure()
plt.imshow(accumulator, extent=[np.rad2deg(theta_range.min()), np.rad2deg(theta_range.max()), rhos.min(), rhos.max()], aspect='auto')
plt.xlabel('Theta (degrees)')
plt.ylabel('Rho')
plt.title('Hough Space with Peaks')
plt.colorbar(label='Votes')
plt.show()



# Plot the detected lines on the noisy data
plt.figure()
plt.scatter(x_noise, y_noise, color='b', label='Noisy data')

for rho_idx, theta_idx in peaks:
    rho = rhos[rho_idx]
    theta = theta_range[theta_idx]
    
    # Convert rho, theta back to Cartesian line (y = mx + c form)
    if np.sin(theta) != 0:
        slope = -np.cos(theta) / np.sin(theta)
        intercept = rho / np.sin(theta)
        x_vals = np.array([x_noise.min(), x_noise.max()])
        y_vals = slope * x_vals + intercept
        plt.plot(x_vals, y_vals, label=f'Line: θ={np.rad2deg(theta):.1f}°', color='r')
    
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.show()