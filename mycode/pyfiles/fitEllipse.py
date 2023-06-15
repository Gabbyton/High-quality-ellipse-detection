
import numpy as np

def drawEllipses(ellipses_para, im):
    x = ellipses_para[0, :].flatten()  # change session 1
    y = ellipses_para[1, :].flatten()

    # Build design matrix
    D = np.column_stack((x * x, x * y, y * y, x, y, np.ones_like(x)))

    # Build scatter matrix
    S = np.dot(D.T, D)

    # Build 6x6 constraint matrix
    C = np.zeros((6, 6))
    C[0, 2] = 2
    C[1, 1] = -1
    C[2, 0] = 2

    # Solve eigensystem
    geval, gevec = np.linalg.eig(S, C)

    # Find the positive eigenvalue
    I = np.where((np.real(np.diag(geval)) > -1e-8) & (~np.isinf(np.diag(geval))))[0]
    if len(I) == 0:
        info = 0
        ellipse = []
        return

    # Extract eigenvector corresponding to negative eigenvalue
    A = np.real(gevec[:, I])

    # unnormalize
    par = A.flatten()  # change session 2

    # Convert to geometric radii and centers
    thetarad = 0.5 * np.arctan2(par[1], par[0] - par[2])
    cost = np.cos(thetarad)
    sint = np.sin(thetarad)
    sin_squared = sint * sint
    cos_squared = cost * cost
    cos_sin = sint * cost
    Ao = par[5]
    Au = par[3] * cost + par[4] * sint
    Av = -par[3] * sint + par[4] * cost
    Auu = par[0] * cos_squared + par[2] * sin_squared + par[1] * cos_sin
    Avv = par[0] * sin_squared + par[2] * cos_squared - par[1] * cos_sin
    tuCentre = -Au / (2 * Auu)
    tvCentre = -Av / (2 * Avv)
    wCentre = Ao - Auu * tuCentre * tuCentre - Avv * tvCentre * tvCentre
    uCentre = tuCentre * cost - tvCentre * sint
    vCentre = tuCentre * sint + tvCentre * cost
    Ru = np.sqrt(np.abs(-wCentre / Auu))
    Rv = np.sqrt(np.abs(-wCentre / Avv))

    # Code by Alan Lu
    # delta   = 4.*par(1).*par(3)-par(2).*par(2);
    # uCentre = (par(2).*par(5)-2.*par(3).*par(4))./delta;
    # vCentre = (par(2).*par(4)-2.*par(1).*par(5))./delta;
    # temp1 = 2.*(par(1).*uCentre.*uCentre+par(3).*vCentre.*vCentre+par(2).*uCentre.*vCentre-par(6));
    # temp2 = par(1)+par(3);
    # temp3 = sqrt((par(1)-par(3)).*(par(1)+par(3))+par(2).*par(2));
    # Ru =
