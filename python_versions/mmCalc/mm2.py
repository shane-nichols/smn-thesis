#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
More efficient partial wave calculation for a single layer. Accepts arrays of
azimuth (Theta), AOI (Phi), and the tensors epilson and alpha. This was written
to efficiently compute MM's of many epsilon and alphas along different
directions and wavelengths for training maching learning algorithms. 
Both ambient media are assumed to have a refractive index of 1.
This calculation uses no for-loops.
@author: Shane Nichols
"""
import numpy as np


def partialWaveFast(e, a, eul, Phi, Theta, d, wavelengths, bReflect):
    # e: electric permittivity. array of shape (N,3,3) with N=len(wavelengths)
    # e: optical activity. array of shape (N,3,3) with N=len(wavelengths)
    #   For one wavelength, e and a can have shape (3,3)
    # eul: array or list of ZXZ Euler rotations in RADIANS to apply to e and a.
    #   These angles correspond to the orientation of the crystal when theta=0
    # Phi: scalar or 1D array of incident angles in RADIANS
    # Theta: scalar or 1D array of azimuthal angles, in RADIANS
    # d: thickness of crystal in nanometers
    # wavelengths: scalar or 1D array of wavelengths in nanometers.
    # bReflect: bool. True for reflection.
    Nwl = np.size(wavelengths)
    Nph = np.size(Phi)
    Naz = np.size(Theta)
    k0 = 2j * np.pi * d / wavelengths
    Ainv = np.array([[0.5, 0.5, 0, 0],
                     [0, 0, 0.5, -0.5j],
                     [0, 0, 0.5, 0.5j],
                     [0.5, -0.5, 0, 0]])
    A = np.array([[1, 0, 0, 1],
                  [1, 0, 0, -1],
                  [0, 1, 1, 0],
                  [0, 1j, -1j, 0]])
    R = R_ZXZ(eul[0] + Theta, eul[1], eul[2]).reshape(Naz, 1, 3, 3)
    if a is 0:
        a = np.broadcast_to(0, [Nph, Naz, Nwl, 3, 3])
    else:
        a = R.swapaxes(-1, -2) @ (1j * a.reshape(1, Nwl, 3, 3)) @ R
        a = np.broadcast_to(a, (Nph, Naz, Nwl, 3, 3))
    e = R.swapaxes(-1, -2) @ e.reshape(1, Nwl, 3, 3) @ R
    e = np.broadcast_to(e, (Nph, Naz, Nwl, 3, 3))
    m = np.broadcast_to(np.eye(3), (Nph, Naz, Nwl, 3, 3))
    kx = np.broadcast_to(np.sin(Phi).reshape(Nph, 1, 1), (Nph, Naz, Nwl))
    idx = [2, 3, 4, 0, 1]  # reindexing list
    S11 = np.array([[a[..., 0, 1], m[..., 1, 1]],
                    [e[..., 0, 0], -a[..., 0, 1]]]).transpose(idx)
    S12 = np.array([[a[..., 1, 1], -m[..., 1, 0]],
                    [e[..., 0, 1], a[..., 0, 0]]]).transpose(idx)
    S21 = np.array([[-a[..., 0, 0], -m[..., 0, 1]],
                    [e[..., 1, 0], -a[..., 1, 1]]]).transpose(idx)
    S22 = np.array([[-a[..., 1, 0], m[..., 0, 0]],
                    [e[..., 1, 1], a[..., 1, 0]]]).transpose(idx)
    S33 = np.array([[-a[..., 2, 2], -m[..., 2, 2]],
                    [-e[..., 2, 2], a[..., 2, 2]]])
    S33 = (S33 / (e[..., 2, 2] * m[..., 2, 2] +
                  a[..., 2, 2] ** 2)).transpose(idx)
    S31 = np.array([[-a[..., 0, 2], -m[..., 2, 1]],
                    [-e[..., 2, 0], a[..., 2, 1] - kx]]).transpose(idx)
    S32 = np.array([[-a[..., 1, 2] + kx, m[..., 2, 0]],
                    [-e[..., 2, 1], -a[..., 2, 0]]]).transpose(idx)
    S13 = np.array([[a[..., 2, 1] + kx, m[..., 1, 2]],
                    [e[..., 0, 2], -a[..., 0, 2]]]).transpose(idx)
    S23 = np.array([[-a[..., 2, 0], -m[..., 0, 2]],
                    [e[..., 1, 2], -a[..., 1, 2] - kx]]).transpose(idx)
    t1 = S33 @ S31
    t2 = S33 @ S32
    delta = np.concatenate(
            (np.concatenate(
                    (S11 - S13 @ t1, S12 - S13 @ t2), axis=4),
             np.concatenate(
                    (S21 - S23 @ t1, S22 - S23 @ t2), axis=4)), axis=3)
    K, Psi1 = np.linalg.eig(delta)
    idx = (np.argsort(K)[..., [3, 2, 0, 1]]).transpose(3, 0, 1, 2)
    Psi1 = Psi1.transpose(4, 3, 0, 1, 2)
    Psi1 = Psi1[(idx[:, None, ...],) + tuple(np.ogrid[:4, :Nph, :Naz, :Nwl])]
    Psi1 = Psi1.transpose(2, 3, 4, 1, 0)
    K = K.transpose(3, 0, 1, 2)
    K = K[(idx, ) + tuple(np.ogrid[:Nph, :Naz, :Nwl])]
    K = K.transpose(1, 2, 3, 0)
    t0 = np.broadcast_to(0, Nph)
    t1 = np.broadcast_to(1, Nph)
    t2 = np.cos(Phi)
    Psi0 = np.array(
            [[t2, t0, -t2, t0],
             [t1, t0, t1, t0],
             [t0, t1, t0, t1],
             [t0, t2, t0, -t2]]). \
        reshape(4, 4, Nph, 1, 1).transpose(2, 3, 4, 0, 1)
    # Parial Wave calculation
    P1 = np.apply_along_axis(
            np.diag, -1, np.exp(-K[..., 0:2] * np.expand_dims(k0, 1)))
    P2 = np.apply_along_axis(
            np.diag, -1, np.exp(K[..., 2:4] * np.expand_dims(k0, 1)))
    Int01 = np.linalg.inv(Psi0) @ Psi1
    Int10 = np.linalg.inv(Psi1) @ Psi0
    T01 = np.linalg.inv(Int01[..., [[0], [1]], [0, 1]])
    T12 = np.linalg.inv(Int10[..., [[0], [1]], [0, 1]])
    R12 = Int10[..., [[2], [3]], [0, 1]] @ T12
    T10 = np.linalg.inv(Int10[..., [[2], [3]], [2, 3]])
    R10 = Int10[..., [[0], [1]], [2, 3]] @ T10
    P1 = bigKron(P1)
    P2 = bigKron(P2)
    T12 = bigKron(T12)
    R12 = bigKron(R12)
    T10 = bigKron(T10)
    R10 = bigKron(R10)

    if bReflect:
        R01 = Int01[..., [[2], [3]], [0, 1]] @ T01
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


def bigKron(a):
    ''' Perfors kron(m, conj(m)) for array of 2x2 matrices.
        a is an N-dimensional array of size N,M,O...2,2
        out is an array of size N,M,O... 4,4'''
    out = np.empty([*a.shape[0:-2], 4, 4], dtype=np.cfloat)
    out[..., 0, 0] = a[..., 0, 0] * np.conj(a[..., 0, 0])
    out[..., 0, 1] = a[..., 0, 0] * np.conj(a[..., 0, 1])
    out[..., 0, 2] = np.conj(out[..., 0, 1])
    out[..., 0, 3] = a[..., 0, 1] * np.conj(a[..., 0, 1])
    out[..., 1, 0] = a[..., 0, 0] * np.conj(a[..., 1, 0])
    out[..., 1, 1] = a[..., 0, 0] * np.conj(a[..., 1, 1])
    out[..., 1, 2] = a[..., 0, 1] * np.conj(a[..., 1, 0])
    out[..., 1, 3] = a[..., 0, 1] * np.conj(a[..., 1, 1])
    out[..., 2, 0] = np.conj(out[..., 1, 0])
    out[..., 2, 1] = np.conj(out[..., 1, 2])
    out[..., 2, 2] = np.conj(out[..., 1, 1])
    out[..., 2, 3] = np.conj(out[..., 1, 3])
    out[..., 3, 0] = a[..., 1, 0] * np.conj(a[..., 1, 0])
    out[..., 3, 1] = a[..., 1, 0] * np.conj(a[..., 1, 1])
    out[..., 3, 2] = np.conj(out[..., 3, 1])
    out[..., 3, 3] = a[..., 1, 1] * np.conj(a[..., 1, 1])
    return out


def R_ZXZ(t1, t2, t3):
    # t1 may be a vector of length N, in which case R has shape N,3,3
    # R_ZXZ
    #    OUT1 = R_ZXZ(T1,T2,T3)
    #    Active ZXZ Euler rotation matrix
    t5 = np.broadcast_to(np.cos(t3), t1.shape)
    t6 = np.sin(t1)
    t7 = np.cos(t1)
    t8 = np.broadcast_to(np.cos(t2), t1.shape)
    t9 = np.broadcast_to(np.sin(t3), t1.shape)
    t10 = np.broadcast_to(np.sin(t2), t1.shape)
    R = np.array([t5 * t7-t6 * t8 * t9, -t7 * t9-t5 * t6 * t8,
                  t6 * t10, t5 * t6 + t7 * t8 * t9,
                 -t6 * t9 + t5 * t7 * t8, -t7 * t10, t9 * t10,
                  t5 * t10, t8]).T.reshape([-1, 3, 3])
    return R
