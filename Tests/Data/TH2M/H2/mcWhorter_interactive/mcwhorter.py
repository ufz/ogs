# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt
import numpy as np

eps2 = 10**-16


class BrooksCorey:
    def __init__(self, pd=5000, lambda_=3, Snr=0.0, Swr=0.0):
        self.pd = pd
        self.lambda_ = lambda_
        self.Snr = Snr
        self.Swr = Swr

    def Sw2Swe(self, Sw):
        return np.clip((Sw - self.Swr) / (1 - self.Snr - self.Swr), eps2, 1 - eps2)

    def Sw2Swe_(func):
        def wrapper(self, Sw):
            Swe = self.Sw2Swe(Sw)
            return func(self, Swe)

        return wrapper

    @Sw2Swe_
    def pc(self, Sw):
        return self.pd * Sw ** (-1 / self.lambda_)

    @Sw2Swe_
    def dPcdSw(self, Sw):
        scaling_factor = 1 / (1 - self.Snr - self.Swr)
        return -self.pd / self.lambda_ * scaling_factor * Sw ** (-1 / self.lambda_ - 1)

    # Relative permeabilities
    @Sw2Swe_
    def kn(self, Sw):
        return (1 - Sw) ** 2 * (1 - Sw ** ((2 + self.lambda_) / self.lambda_))

    @Sw2Swe_
    def kw(self, Sw):
        return Sw ** ((2 + 3 * self.lambda_) / self.lambda_)

    @Sw2Swe_
    def Sw(self, pc):
        return (self.pd / pc) ** self.lambda_

    def plot(self):
        plot_model(self.pc, self.kn, self.kw, Smin=self.Swr, Smax=1 - self.Snr)


class VanGenuchten:
    def __init__(self, alpha=0.0001, n=10, Snr=0.0, Swr=0.0):
        self.alpha = alpha
        self.n = n
        self.m = 1 - 1 / n
        self.Snr = Snr
        self.Swr = Swr

    def Sw2Swe(self, Sw):
        eps2 = 10**-16
        return np.clip((Sw - self.Swr) / (1 - self.Snr - self.Swr), eps2, 1 - eps2)

    def Sw2Swe_(func):
        def wrapper(self, Sw):
            Swe = self.Sw2Swe(Sw)
            return func(self, Swe)

        return wrapper

    @Sw2Swe_
    def pc(self, Sw):
        return (1 / self.alpha) * (Sw ** (-1 / self.m) - 1) ** (1 / self.n)

    @Sw2Swe_
    def dPcdSw(self, Sw):
        scaling_factor = 1 / (1 - self.Snr - self.Swr)
        return (
            -scaling_factor
            / (self.alpha * Sw * self.n * self.m)
            * Sw ** (-1 / self.m)
            * (Sw ** (-1 / self.m) - 1) ** (1 / self.n - 1)
        )

    # Relative permeabilities
    @Sw2Swe_
    def kn(self, Sw):
        return (1 - Sw) ** (1 / 3) * (1 - Sw ** (1 / self.m)) ** (2 * self.m)

    @Sw2Swe_
    def kw(self, Sw):
        return np.sqrt(Sw) * (1 - (1 - Sw ** (1 / self.m)) ** self.m) ** 2

    def plot(self):
        plot_model(self.pc, self.kn, self.kw, Smin=self.Swr, Smax=1 - self.Snr)


def plot_model(pc, kn, kw, Smin=0, Smax=1, nel=100):
    Sw = np.linspace(Smin + eps2, Smax - eps2, nel)

    plt.subplot(1, 2, 1)
    plt.semilogy(Sw, pc(Sw))
    plt.xlabel(r"$S_\mathrm{w}$")
    plt.ylabel(r"$P_\mathrm{c}$")

    ax2 = plt.subplot(1, 2, 2)
    plt.plot(Sw, kn(Sw))
    plt.plot(Sw, kw(Sw))
    plt.xlabel(r"$S_\mathrm{w}$")
    plt.ylabel(r"$k_\mathrm{n}$, $k_\mathrm{w}$")

    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()

    plt.show()


class McWhorter:
    def __init__(
        self,
        model,
        phi=0.5,
        K=1.0e-10,
        muw=0.001,
        mun=0.005,
        S0=0.9,
        Si=0.0,
        t=1000,
        nel=1000,
        max_iter=100000,
        eps=1.0e-14,
    ):
        if abs(S0 - Si) < 1e-8:
            self.A = 0
            self.x = np.linspace(0, 1, nel)
            self.Sw = S0 * np.ones(nel)
        else:
            # Function f and D (McWhorter and Sunada)
            def f(Sw):
                return 1 / (1 + (model.kn(Sw) * muw) / (model.kw(Sw) * mun))

            def D(Sw):
                return -(K * model.kn(Sw) * f(Sw)) / mun * model.dPcdSw(Sw)

            self.f = f
            self.D = D

            def integral(x):
                return np.sum(x)

            Sw = np.linspace(Si + eps, S0 - eps, nel)
            dSw = Sw[1] - Sw[0]

            # Iterative computation of F
            F = np.ones(nel)
            F[0] = eps2
            Falt = F.copy()

            Di = D(Sw)

            for _it in range(max_iter):
                a = integral((Sw[1:] - Si) * Di[1:] / F[1:])

                for i in range(1, len(Sw) - 1):
                    b = integral((Sw[i:] - Sw[i]) * Di[i:] / F[i:])
                    F[i] = 1 - b / a

                if np.linalg.norm(F - Falt) < eps:
                    break

                Falt = F.copy()
                self.A = np.sqrt(phi / 2 * a * dSw)

            dFdS = np.zeros_like(F)
            for i in range(len(Sw)):
                dFdS[i] = integral(Di[i:] / F[i:]) / a

            x = 2 * self.A / phi * dFdS * np.sqrt(t)

            self.x = np.flip(x[1:])  # x values must be increasing for np.interp
            self.Sw = np.flip(Sw[1:])

    def get_solution(self):
        return [self.x, self.Sw]

    def plot_solution(self):
        plt.plot(self.x, self.Sw, label=f"Analytical solution A={self.A:.4e}")
        plt.xlabel(r"$x$ (m)")
        plt.ylabel(r"$S_\mathrm{w}$ (-)")
        plt.legend()
        # plt.show()

    def plot_D_f(self):
        Sw = self.Sw

        plt.subplot(1, 2, 1)
        plt.plot(Sw, self.D(Sw))
        plt.xlabel(r"$S_\mathrm{w}$")
        plt.ylabel(r"$D$")

        ax2 = plt.subplot(1, 2, 2)
        plt.plot(Sw, self.f(Sw))
        plt.xlabel(r"$S_{w}$")
        plt.ylabel(r"$f$")

        ax2.yaxis.set_label_position("right")
        ax2.yaxis.tick_right()

        plt.show()
