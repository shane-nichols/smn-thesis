#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Matlab function mmPartialWave and its dependencies converted to Python

@author: Shane Nichols
"""
import numpy as np


def mmPartialWave(layerArray, wavelengths, aoi, bReflect, bNorm):
    aoi = aoi * np.pi / 180
    Nlayers = len(layerArray) - 1
    g = [i for i, elem in enumerate(layerArray) if elem[3] is True][0]
    if layerArray[0][0] == 'air':
        t1 = np.cos(aoi)
        Psi0 = np.array([[t1, 0, -t1, 0],
                         [1, 0, 1, 0],
                         [0, 1, 0, 1],
                         [0, t1, 0, -t1]])
        n0 = np.ones(len(wavelengths))
    else:
        Psi0, n0 = psiAmbient(layerArray[0], aoi, wavelengths)

    kx = n0 * np.sin(aoi)  # x-component of the wavevector

    if layerArray[Nlayers][0] == 'air':
        t1 = np.cos(aoi)
        Psi2 = np.array([[t1, 0, -t1, 0],
                         [1, 0, 1, 0],
                         [0, 1, 0, 1],
                         [0, t1, 0, -t1]])
    elif layerArray[Nlayers][4] == 0:
        Psi2 = psiAniso(layerArray[Nlayers], wavelengths, kx)
    elif layerArray[Nlayers][1] == np.inf:
        Psi2 = psiIso(layerArray[Nlayers], wavelengths, kx)
    else:
        Psi2 = psiAmbient(layerArray[Nlayers], aoi, wavelengths)

    # if there are thin layers before thick one, compute layer matrices
    if g > 1:
        for m in np.arange(1, g, 1):
            # change sign of d to invert layer matrix
            layerArray[m][1] = -layerArray[m][1]
            if layerArray[m][4] == 0:  # check if layer is anisotropic
                Psi0 = layerBerreman(layerArray[m], wavelengths, kx) @ Psi0
            else:
                Psi0 = layerBerremanIso(layerArray[m], wavelengths, kx) @ Psi0

    # if there are thin layers after thick one, compute layer matrices
    if Nlayers - g > 1:
        for m in np.arange(Nlayers - 1, g, -1):
            if layerArray[m][4] == 0:   # check if layer is anisotropic
                Psi2 = layerBerreman(layerArray[m], wavelengths, kx) @ Psi2
            else:
                Psi2 = layerBerremanIso(layerArray[m], wavelengths, kx) @ Psi2

    M = partialWave(Psi0, Psi2, layerArray[g], wavelengths, kx, bReflect)

    if bNorm:  # normalize by M11 is bNorm is set
        M = M.transpose([1, 2, 0])
        M = M / M[0, 0, :]
        M = M.transpose([2, 0, 1])
    return M


def psiAmbient(layer, aoi, wavelengths):
    # aoi : free-space angle of incidence in radians
    # layer : layer containing the material name
    n = materialLibIso(layer[0], wavelengths)
    aoi = np.cos(aoi)
    Psi = np.zeros([len(wavelengths), 4, 4])
    for j in range(len(wavelengths)):
        Psi[j] = np.array([[aoi, 0, -aoi, 0],
                           [n[j], 0, n[j], 0],
                           [0, 1, 0, 1],
                           [0, n[j] * aoi, 0, -n[j] * aoi]])
    return Psi, n


def psiAniso(layer, wavelengths, kx):
    epsilon, alpha, mu = materialLib(layer[0], wavelengths)
    eul = layer[2] * np.pi / 180
    Nwl = len(wavelengths)
    R = R_ZXZ(eul[0], eul[1], eul[2])
    if alpha == 0:
        alpha = np.zeros([Nwl, 3, 3])
    else:
        alpha = R.T @ (1j * alpha) @ R
    epsilon = R.T @ epsilon @ R
    mu = R.T @ mu @ R
    delta = deltaMatAr(epsilon, alpha, mu, kx)
    K, Psi = np.linalg.eig(delta)
    idx = np.argsort(K)[:, [3, 2, 0, 1]]
    for i in range(Nwl):
        K[i] = K[i, idx[i, :]]
        Psi[i] = Psi[i, :, idx[i, :]].T
    return Psi


def psiIso(layer, wavelengths, kx):
    epsilon = materialLibIso(layer[0], wavelengths) ** 2
    psi = np.zeros([len(wavelengths), 4, 4], dtype=np.cfloat)
    q = np.sqrt(epsilon - kx ** 2)
    t1 = q / epsilon
    t2 = 1 / q
    t3 = np.ones([1, len(wavelengths)])
    psi[:, 0, 0] = t1
    psi[:, 0, 3] = -t2
    psi[:, 1, 0] = t3
    psi[:, 1, 3] = t3
    psi[:, 2, 1] = t2
    psi[:, 2, 2] = -t1
    psi[:, 3, 1] = t3
    psi[:, 3, 2] = t3
    return psi


def layerBerreman(layer, wavelengths, kx):
    # calculates the INVERSE Berreman layer matrix, i.e., the quantity L(-d) in
    # J. Opt. Soc. Am. A, 32, 2015, 2049-2057.

    epsilon, alpha, mu = materialLib(layer[0], wavelengths)
    d = layer[1]
    k0 = 2j * np.pi * d / wavelengths
    eul = layer[2] * np.pi/180
    Nwl = len(wavelengths)
    R = R_ZXZ(eul[0], eul[1], eul[2])
    if alpha == 0:
        alpha = np.zeros([Nwl, 3, 3])
    else:
        alpha = R.T @ (1j * alpha) @ R
    epsilon = R.T @ epsilon @ R
    mu = R.T @ mu @ R
    delta = deltaMatAr(epsilon, alpha, mu, kx)
    K, Psi = np.linalg.eig(delta)
    layerMatrix = Psi @ \
        setDiag(np.exp(K * np.expand_dims(k0, 1))) @ np.linalg.inv(Psi)
    return layerMatrix


def layerBerremanIso(layer, wavelengths, kx):
    # calculates the INVERSE Berreman layer matrix for isotropic media.
    # i.e., the quantity L(-d) in J. Opt. Soc. Am. A, 32, 2015, 2049-2057.

    epsilon = materialLibIso(layer[0], wavelengths) ** 2
    d = layer[1]
    q = np.sqrt(epsilon - kx ** 2)
    arg = (2 * np.pi * d) * (q / wavelengths)
    C = np.cos(arg)
    S = 1j * np.sin(arg)
    qS = q * S
    Sq = S / q
    layerMatrix = np.zeros([len(wavelengths), 4, 4], dtype=np.cfloat)
    layerMatrix[:, 0, 0] = C
    layerMatrix[:, 0, 1] = qS / epsilon
    layerMatrix[:, 1, 0] = epsilon * Sq
    layerMatrix[:, 1, 1] = C
    layerMatrix[:, 2, 2] = C
    layerMatrix[:, 2, 3] = Sq
    layerMatrix[:, 3, 2] = qS
    layerMatrix[:, 3, 3] = C
    return layerMatrix


def partialWave(Psi0, Psi2, layer, wavelengths, kx, bReflect):

    '''epsilon and alpha: 3x3 constitutive tensors in the standard setting
       eul: array of ZXZ passive euler angles to rotate the tensors, in radians
       AOI: incident angle, in radians '''
    epsilon, alpha, mu = materialLib(layer[0], wavelengths)
    d = layer[1]
    k0 = 2j * np.pi * d / wavelengths
    eul = layer[2] * np.pi / 180
    Nwl = len(wavelengths)
    Ainv = np.array([[0.5, 0.5, 0, 0],
                     [0, 0, 0.5, -0.5j],
                     [0, 0, 0.5, 0.5j],
                     [0.5, -0.5, 0, 0]])
    A = np.array([[1, 0, 0, 1],
                  [1, 0, 0, -1],
                  [0, 1, 1, 0],
                  [0, 1j, -1j, 0]])
    R = R_ZXZ(eul[0], eul[1], eul[2])
    if alpha is 0:
        alpha = np.zeros([Nwl, 3, 3])
    else:
        alpha = R.T @ (1j * alpha) @ R
    epsilon = R.T @ epsilon @ R
    mu = R.T @ mu @ R
    delta = deltaMatAr(epsilon, alpha, mu, kx)
    K, Psi1 = np.linalg.eig(delta)
    idx = np.argsort(K)[:, [3, 2, 0, 1]]
    for i in range(Nwl):
        K[i] = K[i, idx[i, :]]
        Psi1[i] = Psi1[i, :, idx[i, :]].T
    # Parial Wave calculation
    Psi1inv = np.linalg.inv(Psi1)
    P1 = setDiag(np.exp(-K[:, 0:2] * np.expand_dims(k0, 1)))
    P2 = setDiag(np.exp(K[:, 2:4] * np.expand_dims(k0, 1)))
    Int01 = np.linalg.inv(Psi0) @ Psi1
    Int12 = Psi1inv @ Psi2
    Int10 = Psi1inv @ Psi0
    T01 = np.linalg.inv(Int01[:, [[0], [1]], [0, 1]])
    T12 = np.linalg.inv(Int12[:, [[0], [1]], [0, 1]])
    R12 = Int12[:, [[2], [3]], [0, 1]] @ T12
    T10 = np.linalg.inv(Int10[:, [[2], [3]], [2, 3]])
    R10 = Int10[:, [[0], [1]], [2, 3]] @ T10
    P1 = bigKron(P1)
    P2 = bigKron(P2)
    T12 = bigKron(T12)
    R12 = bigKron(R12)
    T10 = bigKron(T10)
    R10 = bigKron(R10)

    if bReflect:
        R01 = Int01[:, [[2], [3]], [0, 1]] @ T01
        T01 = bigKron(T01)
        R01 = bigKron(R01)
        P2 = P2 @ R12 @ P1  # redefine P2 for in-place
        MM = R01 + T10 @ P2 @ np.linalg.inv(np.eye(4) - R10 @ P2) @ T01
        MM = np.real(A @ MM @ Ainv)
    else:
        T01 = bigKron(T01)
        MM = T12 @ P1 @ np.linalg.inv(np.eye(4) - R10 @ P2 @ R12 @ P1) @ T01
        MM = np.real(A @ MM @ Ainv)
    return MM


def deltaMat(M, kx):  # for element-by-element construction using 6x6 matrix M
    # no longer used because it's crazy slow
    S11 = np.array(-M[[[4], [0]], [0, 4]])
    S12 = np.array([[-M[4, 1],  M[4, 3]],
                    [-M[0, 1],  M[0, 3]]])
    S21 = np.array([[M[3, 0],  M[3, 4]],
                    [-M[1, 0], -M[1, 4]]])
    S22 = np.array([[M[3, 1], -M[3, 3]],
                    [-M[1, 1], M[1, 3]]])
    S33 = np.array([[-M[2, 5],  M[5, 5]],
                    [M[2, 2],  -M[5, 2]]]) \
        / (M[2, 2] * M[5, 5] - M[2, 5] * M[5, 2])
    S31 = np.array([[M[5, 0],  M[5, 4]],
                    [M[2, 0],  M[2, 4] + kx]])
    S32 = np.array([[M[5, 1] - kx, -M[5, 3]],
                    [M[2, 1],         -M[2, 3]]])
    S13 = -np.array([[M[4, 2] + kx, M[4, 5]],
                     [M[0, 2],         M[0, 5]]])
    S23 = np.array([[M[3, 2],   M[3, 5]],
                    [-M[1, 2], -M[1, 5] + kx]])
    t1 = S33 @ S31
    t2 = S33 @ S32
    delta = -np.bmat([[S11 - S13 @ t1, S12 - S13 @ t2],
                      [S21 - S23 @ t1, S22 - S23 @ t2]])
    return delta


def deltaMatAr(e, a, m, kx):
    # build Delta all at once using the 3x3 constitutive tensors
    # e is epsilon
    # a is alpha
    # m is mu
    S11 = np.array([[a[:, 0, 1], m[:, 1, 1]],
                    [e[:, 0, 0], -a[:, 0, 1]]]).transpose([2, 0, 1])
    S12 = np.array([[a[:, 1, 1], -m[:, 1, 0]],
                    [e[:, 0, 1], a[:, 0, 0]]]).transpose([2, 0, 1])
    S21 = np.array([[-a[:, 0, 0], -m[:, 0, 1]],
                    [e[:, 1, 0], -a[:, 1, 1]]]).transpose([2, 0, 1])
    S22 = np.array([[-a[:, 1, 0], m[:, 0, 0]],
                    [e[:, 1, 1], a[:, 1, 0]]]).transpose([2, 0, 1])
    S33 = np.array([[-a[:, 2, 2], -m[:, 2, 2]],
                    [-e[:, 2, 2], a[:, 2, 2]]])
    S33 = (S33 /
           (e[:, 2, 2] * m[:, 2, 2] + a[:, 2, 2] ** 2)).transpose([2, 0, 1])
    S31 = np.array([[-a[:, 0, 2], -m[:, 2, 1]],
                    [-e[:, 2, 0], a[:, 2, 1] - kx]]).transpose([2, 0, 1])
    S32 = np.array([[-a[:, 1, 2] + kx, m[:, 2, 0]],
                    [-e[:, 2, 1], -a[:, 2, 0]]]).transpose([2, 0, 1])
    S13 = np.array([[a[:, 2, 1] + kx, m[:, 1, 2]],
                    [e[:, 0, 2], -a[:, 0, 2]]]).transpose([2, 0, 1])
    S23 = np.array([[-a[:, 2, 0], -m[:, 0, 2]],
                    [e[:, 1, 2], -a[:, 1, 2] - kx]]).transpose([2, 0, 1])
    t1 = S33 @ S31
    t2 = S33 @ S32
    delta = np.concatenate(
                (np.concatenate(
                        (S11 - S13 @ t1, S12 - S13 @ t2), axis=2),
                 np.concatenate(
                        (S21 - S23 @ t1, S22 - S23 @ t2), axis=2)), axis=1)
    return delta


def setDiag(diag):
    ''' sets the diagonals of an array of matrices.
    diag is a 2D array where each column is the diagonal elements of a matrix.
    out is a 3D array of diagonal matrices. Although I use a loop, I loop over
    the diagonal elements, not the matrix index, so it is fast'''
    out = np.zeros([diag.shape[0], diag.shape[1], diag.shape[1]],
                   dtype=np.cfloat)
    for i in range(diag.shape[1]):
        out[:, i, i] = diag[:, i]
    return out


def bigKron(a):
    ''' Perfors kron(m, conj(m)) for array of 2x2 matrices.
        a is an N-dimensional array of size N,M,O...2,2
        out is an array of size N,M,O... 4,4'''
    sz = a.shape
    a = a.reshape([-1, 2, 2])
    out = np.zeros([a.shape[0], 4, 4], dtype=np.cfloat)
    out[:, 0, 0] = a[:, 0, 0] * np.conj(a[:, 0, 0])
    out[:, 0, 1] = a[:, 0, 0] * np.conj(a[:, 0, 1])
    out[:, 0, 2] = np.conj(out[:, 0, 1])
    out[:, 0, 3] = a[:, 0, 1] * np.conj(a[:, 0, 1])
    out[:, 1, 0] = a[:, 0, 0] * np.conj(a[:, 1, 0])
    out[:, 1, 1] = a[:, 0, 0] * np.conj(a[:, 1, 1])
    out[:, 1, 2] = a[:, 0, 1] * np.conj(a[:, 1, 0])
    out[:, 1, 3] = a[:, 0, 1] * np.conj(a[:, 1, 1])
    out[:, 2, 0] = np.conj(out[:, 1, 0])
    out[:, 2, 1] = np.conj(out[:, 1, 2])
    out[:, 2, 2] = np.conj(out[:, 1, 1])
    out[:, 2, 3] = np.conj(out[:, 1, 3])
    out[:, 3, 0] = a[:, 1, 0] * np.conj(a[:, 1, 0])
    out[:, 3, 1] = a[:, 1, 0] * np.conj(a[:, 1, 1])
    out[:, 3, 2] = np.conj(out[:, 3, 1])
    out[:, 3, 3] = a[:, 1, 1] * np.conj(a[:, 1, 1])
    out = out.reshape(np.concatenate((sz[0:len(sz)-2], [4, 4])))
    return out


def R_ZXZ(t1, t2, t3):
    # R_ZXZ
    #    OUT1 = R_ZXZ(T1,T2,T3)
    #    Active ZXZ Euler rotation matrix
    t5 = np.cos(t3)
    t6 = np.sin(t1)
    t7 = np.cos(t1)
    t8 = np.cos(t2)
    t9 = np.sin(t3)
    t10 = np.sin(t2)
    R = np.reshape(
            np.array([t5 * t7-t6 * t8 * t9, -t7 * t9-t5 * t6 * t8,
                      t6 * t10, t5 * t6 + t7 * t8 * t9,
                      -t6 * t9 + t5 * t7 * t8, -t7 * t10, t9 * t10,
                      t5 * t10, t8]), [3, 3]).T
    return R


def materialLib(material, wavelengths):

    materialDict = {
            '+quartz': optFun_RHquartz(wavelengths),
            '-quartz': optFun_LHquartz(wavelengths),
            'TiO2': optFun_TiO2(wavelengths)
            }
    return materialDict.get(material)


def optFun_RHquartz(wavelengths):
    Nwl = len(wavelengths)
    epsilon = np.zeros([Nwl, 3, 3])
    mu = setDiag(np.ones([Nwl, 3]))
    alpha = np.zeros([Nwl, 3, 3])
    lam2 = (wavelengths / 1000) ** 2
    epsilon[:, 0, 0] = 1.07044083 * lam2 / (lam2 - 0.0100585997) \
        + 1.10202242 * lam2 / (lam2 - 100) + 1.28604141
    epsilon[:, 1, 1] = epsilon[:, 0, 0]
    epsilon[:, 2, 2] = 1.09509924 * lam2 / (lam2 - 0.0102101864) \
        + 1.15662475 * lam2 / (lam2 - 100) + 1.28851804
    lam2 = wavelengths ** 2
    alpha[:, 0, 0] = wavelengths ** 3 * (0.0198) / (lam2 - 93 ** 2) ** 2
    alpha[:, 2, 2] = wavelengths ** 3 * (-0.0408) / (lam2 - 87 ** 2) ** 2
    alpha[:, 1, 1] = alpha[:, 0, 0]
    return epsilon, alpha, mu


def optFun_LHquartz(wavelengths):
    Nwl = len(wavelengths)
    epsilon = np.zeros([Nwl, 3, 3])
    mu = setDiag(np.ones([Nwl, 3]))
    alpha = np.zeros([Nwl, 3, 3])
    lam2 = (wavelengths / 1000) ** 2
    epsilon[:, 0, 0] = 1.07044083 * lam2 / (lam2 - 0.0100585997) \
        + 1.10202242 * lam2 / (lam2 - 100) + 1.28604141
    epsilon[:, 1, 1] = epsilon[:, 0, 0]
    epsilon[:, 2, 2] = 1.09509924 * lam2 / (lam2 - 0.0102101864) \
        + 1.15662475 * lam2 / (lam2 - 100) + 1.28851804
    lam2 = wavelengths ** 2
    alpha[:, 0, 0] = wavelengths ** 3 * (-0.0198) / (lam2 - 93 ** 2) ** 2
    alpha[:, 2, 2] = wavelengths ** 3 * (0.0408) / (lam2 - 87 ** 2) ** 2
    alpha[:, 1, 1] = alpha[:, 0, 0]
    return epsilon, alpha, mu


def optFun_TiO2(wavelengths):
    Nwl = len(wavelengths)
    epsilon = np.zeros([Nwl, 3, 3])
    mu = setDiag(np.ones([Nwl, 3]))
    lam2 = (wavelengths/1000) ** 2
    epsilon[:, 0, 0] = 0.2441 / (lam2 - 0.0803) + 5.913
    epsilon[:, 1, 1] = epsilon[:, 0, 0]
    epsilon[:, 2, 2] = 0.3322 / (lam2 - 0.0843) + 7.197
    return epsilon, 0, mu


def materialLibIso(material, wavelengths):

    materialDict = {
            'air': np.ones(len(wavelengths)),
            'BK7': optFun_BK7(wavelengths)
            }
    return materialDict.get(material)


def optFun_BK7(wavelengths):
    oscA = np.array([1.03961212, 0.231792344, 1.01046945])
    oscE = np.array([6000.69867, 20017.9144, 103560653])
    lam2 = wavelengths ** 2
    n = np.zeros(len(wavelengths))
    for i in range(len(wavelengths)):
        n[i] = np.sqrt(np.sum(lam2[i] * oscA / (lam2[i] - oscE)) + 1)
    return n
