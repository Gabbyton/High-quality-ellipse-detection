import cv2
import numpy as np

def LCS_ellipse(image):
    # Preprocess the image
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    blurred = cv2.GaussianBlur(gray, (5, 5), 0)
    edges = cv2.Canny(blurred, 50, 150)

    # Perform Hough Transform to detect ellipses
    circles = cv2.HoughCircles(edges, cv2.HOUGH_GRADIENT, dp=1, minDist=100,
                               param1=50, param2=30, minRadius=0, maxRadius=0)

    # Convert the detected circles to ellipses
    ellipses = []
    if circles is not None:
        circles = np.round(circles[0, :]).astype(int)
        for (x, y, r) in circles:
            ellipse = cv2.fitEllipse(np.vstack(np.where(edges[y-r:y+r, x-r:x+r] > 0)).T)
            ellipses.append(ellipse)

    return ellipses

# Load an image
image = cv2.imread('your_image.jpg')

# Detect ellipses
ellipses = LCS_ellipse(image)

# Display the detected ellipses
for ellipse in ellipses:
    cv2.ellipse(image, ellipse, (0, 255, 0), 2)

cv2.imshow('Detected Ellipses', image)
cv2.waitKey(0)
cv2.destroyAllWindows()
