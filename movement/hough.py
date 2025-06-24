# -*- coding: utf-8 -*-


import cv2
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

# Generate noisy 2D data
np.random.seed(42)
x = np.linspace(0, 100, 500)
y = 2 * x + 1 + np.random.normal(0, 10, x.shape)

# Create an image (black background)
image_size = 512
img = np.zeros((image_size, image_size), dtype=np.uint8)

# Normalize and add points to the image
x_norm = np.interp(x, (x.min(), x.max()), (0, image_size - 1)).astype(int)
y_norm = np.interp(y, (y.min(), y.max()), (0, image_size - 1)).astype(int)

for xi, yi in zip(x_norm, y_norm):
    img[yi, xi] = 255  # Mark the point with white color (edge)

# Display the image
plt.figure()
plt.imshow(img, cmap='gray')
plt.title("Noisy Points Image")
plt.show()

# Apply Canny Edge Detector to find edges
edges = cv2.Canny(img, 50, 150)

# Display edges
plt.figure()
plt.imshow(edges, cmap='gray')
plt.title("Edges Detected")
plt.show()

# Perform Hough Line Transform
lines = cv2.HoughLines(edges, 2, np.pi / 180, 150)

# Draw the detected lines on the original image
if lines is not None:
    for rho, theta in lines[:, 0]:
        a = np.cos(theta)
        b = np.sin(theta)
        x0 = a * rho
        y0 = b * rho
        x1 = int(x0 + 1000 * (-b))
        y1 = int(y0 + 1000 * (a))
        x2 = int(x0 - 1000 * (-b))
        y2 = int(y0 - 1000 * (a))
        
        cv2.line(img, (x1, y1), (x2, y2), (255, 0, 0), 2)

# Display the image with detected lines
plt.figure()
plt.imshow(img, cmap='gray')
plt.title("Detected Lines with Hough Transform")
plt.show()


