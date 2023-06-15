import numpy as np
import cv2
import matplotlib.pyplot as plt

def drawEllipses(ellipses_para, im):
    if im is not None:
        plt.figure()
        plt.imshow(im, 'border', 'tight', 'initialmagnification', 'fit')
        size_im = im.shape[:2]
        plt.gca().set_aspect('equal')
        plt.gca().set_xlim([0, size_im[1]])
        plt.gca().set_ylim([size_im[0], 0])
        plt.gca().invert_yaxis()
        plt.axis('off')
    else:
        plt.gca().set_aspect('equal')
        plt.axis('off')

    th = np.arange(0, 2 * np.pi, np.pi / 180)
    for i in range(ellipses_para.shape[1]):
        Semi_major = ellipses_para[2, i]
        Semi_minor = ellipses_para[3, i]
        x0 = ellipses_para[0, i]
        y0 = ellipses_para[1, i]
        Phi = ellipses_para[4, i]
        x = x0 + Semi_major * np.cos(Phi) * np.cos(th) - Semi_minor * np.sin(Phi) * np.sin(th)
        y = y0 + Semi_minor * np.cos(Phi) * np.sin(th) + Semi_major * np.sin(Phi) * np.cos(th)
        plt.plot(x, y, 'r', linewidth=2)

    if im is not None:
        plt.show()
