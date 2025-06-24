# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import RANSACRegressor
from sklearn.linear_model import LinearRegression

plt.close('all')

# Step 1: Generate synthetic data
np.random.seed(42)
n_samples = 100
n_outliers = 20

# Generate line data with some noise
X = np.linspace(-10, 10, n_samples).reshape(-1, 1) + 3*np.random.normal(size=n_samples).reshape(-1, 1)
y = 2 * X.squeeze() + 1 + 3*np.random.normal(size=n_samples)

# Add some outliers
X[:n_outliers] = 10 * np.random.normal(size=(n_outliers, 1))
y[:n_outliers] = 20 * np.random.normal(size=n_outliers)

# Step 2: Fit the model using RANSAC
ransac = RANSACRegressor(LinearRegression(), residual_threshold=5.0, random_state=42)
ransac.fit(X, y)
inlier_mask = ransac.inlier_mask_
outlier_mask = ~inlier_mask

# Step 3: Plot the data
line_X = np.arange(X.min(), X.max())[:, np.newaxis]
line_y_ransac = ransac.predict(line_X)

plt.scatter(X[inlier_mask], y[inlier_mask], color="yellowgreen", marker="o", label="Inliers")
plt.scatter(X[outlier_mask], y[outlier_mask], color="red", marker="x", label="Outliers")
plt.plot(line_X, line_y_ransac, color="cornflowerblue", linewidth=2, label="RANSAC fit")
plt.legend(loc="upper left")
plt.xlabel("X")
plt.ylabel("y")
plt.title("RANSAC Line Fitting")
plt.show()