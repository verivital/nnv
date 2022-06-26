# different classes with benchmarks

import jax.numpy as np
from jax.numpy import tanh
from jax.numpy import sin
from jax.numpy import cos
from jax.numpy import exp


def get_model(benchmark, radius=None):
    if benchmark == "bruss":
        return Brusselator(radius)  # Benchmark to run
    elif benchmark == "vdp":
        return VanDerPol(radius)  # Benchmark to run
    elif benchmark == "robot":
        return Robotarm(radius)  # Benchmark to run
    elif benchmark == "dubins":
        return DubinsCar(radius)  # Benchmark to run
    elif benchmark == "ms":
        return MitchellSchaeffer(radius)  # Benchmark to run
    elif benchmark == "cartpole":
        return CartpoleLinear(radius)  # Benchmark to run
    elif benchmark == "quadcopter":
        return Quadcopter(radius)  # Benchmark to run
    elif benchmark == "cartpoleCTRNN":
        return CartpoleCTRNN(radius)  # Benchmark to run
    elif benchmark == "cartpoleLTC":
        return CartpoleLTC(radius)  # Benchmark to run
    elif benchmark == "cartpoleLTC_RK":
        return CartpoleLTC_RK(radius)  # Benchmark to run
    elif benchmark == "ldsCTRNN":
        return LDSwithCTRNN(radius)  # Benchmark to run
    elif benchmark == "pendulumCTRNN":
        return PendulumwithCTRNN(radius)  # Benchmark to run
    elif benchmark == "fpaCTRNN":  # 
        return CTRNN_fpa(radius)
    elif benchmark == "CTRNNosc":
        return CTRNNosc(radius)  # Benchmark to run
    elif benchmark == "spiralL":
        return spiralL(radius)  # Benchmark to run
    elif benchmark == "spiralNL":
        return spiralNL(radius)  # Benchmark to run
    else:
        raise ValueError("Unknown benchmark " + benchmark)


# 2-dimensional brusselator
class Brusselator:
    def __init__(self, radius=None):
        # ============ adapt initial values ===========
        self.cx = (1, 1)
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 0.01
        # ===================================================
        self.cx = np.array(self.cx, dtype=float)
        self.dim = self.cx.size  # dimension of the system

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim, dtype=object)

        # ============ adapt input and system dynamics ===========
        x, y = x

        a = 1
        b = 1.5

        fx = a + x ** 2 * y - (b + 1) * x
        fy = b * x - x ** 2 * y

        system_dynamics = [fx, fy]  # has to be in the same order as the input variables
        # ===================================================

        return np.array(system_dynamics)  # return as numpy array


# 2-dimensional van der pol
class VanDerPol:
    def __init__(self, radius):
        # ============ adapt initial values ===========
        self.cx = (-1, -1)
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 0.5
        # ===================================================
        self.cx = np.array(self.cx, dtype=float)
        self.dim = self.cx.size  # dimension of the system

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim, dtype=object)

        # ============ adapt input and system dynamics ===========
        x, y = x

        fx = y
        fy = (x * x - 1) * y - x

        system_dynamics = [fx, fy]  # has to be in the same order as the input variables
        # ===================================================

        return np.array(system_dynamics)  # return as numpy array


# 4-dimensional Robotarm
class Robotarm:
    def __init__(self, radius):
        # ============ adapt initial values ===========
        self.cx = (1.505, 1.505, 0.005, 0.005)
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 0.05
        # ===================================================
        self.cx = np.array(self.cx, dtype=float)
        self.dim = self.cx.size  # dimension of the system

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim, dtype=object)

        # ============ adapt input and system dynamics ===========
        x1, x2, x3, x4 = x

        m = 1
        l = 3
        kp1 = 2
        kp2 = 1
        kd1 = 2
        kd2 = 1

        fx1 = x3
        fx2 = x4
        fx3 = (-2 * m * x2 * x3 * x4 - kp1 * x1 - kd1 * x3) / (m * x2 * x2 + l / 3) + (
            kp1 * kp1
        ) / (m * x2 * x2 + l / 3)
        fx4 = x2 * x3 * x3 - kp2 * x2 / m - kd2 * x4 / m + kp2 * kp2 / m

        system_dynamics = [
            fx1,
            fx2,
            fx3,
            fx4,
        ]  # has to be in the same order as the input variables
        # ===================================================

        return np.array(system_dynamics)  # return as numpy array


# 3-dimensional dubins car
class DubinsCar:
    def __init__(self, radius):
        # ============ adapt initial values ===========
        self.cx = (0, 0, 0.7854, 0)
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 0.01
        # ===================================================
        self.cx = np.array(self.cx, dtype=float)
        self.dim = self.cx.size  # dimension of the system

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim, dtype=object)

        # ============ adapt input and system dynamics ===========
        x, y, th, tt = x

        v = 1

        fx = v * cos(th)
        fy = v * sin(th)
        fth = x * sin(tt)
        ftt = np.array(1)  # this is needed for lax_numpy to see this as an array

        system_dynamics = [
            fx,
            fy,
            fth,
            ftt,
        ]  # has to be in the same order as the input variables
        # ===================================================

        return np.array(system_dynamics)  # return as numpy array


# 2-dimensional Mitchell  Schaeffer  cardiac-cell
class MitchellSchaeffer:
    def __init__(self, radius):
        # ============ adapt initial values ===========
        self.cx = (0.8, 0.5)
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 0.1
        # ===================================================
        self.cx = np.array(self.cx, dtype=float)
        self.dim = self.cx.size  # dimension of the system

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim, dtype=object)

        # ============ adapt input and system dynamics ===========
        x1, x2 = x

        sig_x1 = 0.5 * (1 + tanh(50 * x1 - 5))

        fx1 = x2 * x1 ** 2 * (1 - x1) / 0.3 - x1 / 6
        fx2 = sig_x1 * (-x2 / 150) + (1 - sig_x1) * (1 - x2) / 20

        system_dynamics = [
            fx1,
            fx2,
        ]  # has to be in the same order as the input variables
        # ===================================================

        return np.array(system_dynamics)  # return as numpy array


# 4-dimensional cartpole with linear stabilizing controller
class CartpoleLinear:
    def __init__(self, radius):
        # ============ adapt initial values ===========
        self.cx = (0, 0, 0.001, 0)  # initial values
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 0.05
        # ===================================================

        self.cx = np.array(self.cx, dtype=float)
        self.dim = self.cx.size  # dimension of the system

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim, dtype=object)

        # ============ adapt input and system dynamics ===========
        dth, dx, th, x = x  # input variables

        M = 1.0
        g = 9.81
        l = 1.0
        m = 0.001

        f = -1.1 * M * g * th - dth

        fdth = (
            1.0
            / (l * (M + m * sin(th) * sin(th)))
            * (
                f * cos(th)
                - m * l * dth * dth * cos(th) * sin(th)
                + (m + M) * g * sin(th)
            )
        )
        fdx = (
            1.0
            / (M + m * sin(th) * sin(th))
            * (f + m * sin(th) * (-l * dth * dth + g * cos(th)))
        )

        fx = dx

        fth = dth

        system_dynamics = [
            fdth,
            fdx,
            fth,
            fx,
        ]  # has to be in the same order as the input variables
        # ===================================================

        return np.array(system_dynamics)  # return as numpy array


# 17-dimensional Quadcopter
class Quadcopter:
    def __init__(self, radius):
        # ============ adapt initial values ===========
        self.cx = (
            -0.995,
            -0.995,
            9.005,
            -0.995,
            -0.995,
            -0.995,
            -0.995,
            -0.995,
            -0.995,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
        )
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 0.005
        # ===================================================
        self.cx = np.array(self.cx, dtype=float)
        self.dim = self.cx.size  # dimension of the system

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim, dtype=object)

        # ============ adapt input and system dynamics ===========
        pn, pe, h, u, v, w, p, q, r, q0, q1, q2, q3, pI, qI, rI, hI = x

        pr = 0
        qr = 0
        rr = 0
        hr = 0

        pn_ = (
            2 * u * (q0 * q0 + q1 * q1 - 0.5)
            - 2 * v * (q0 * q3 - q1 * q2)
            + 2 * w * (q0 * q2 + q1 * q3)
        )
        pe_ = (
            2 * v * (q0 * q0 + q2 * q2 - 0.5)
            + 2 * u * (q0 * q3 + q1 * q2)
            - 2 * w * (q0 * q1 - q2 * q3)
        )
        h_ = (
            2 * w * (q0 * q0 + q3 * q3 - 0.5)
            - 2 * u * (q0 * q2 - q1 * q3)
            + 2 * v * (q0 * q1 + q2 * q3)
        )

        u_ = r * v - q * w - 11.62 * (q0 * q2 - q1 * q3)
        v_ = p * w - r * u + 11.62 * (q0 * q1 + q2 * q3)
        w_ = q * u - p * v + 11.62 * (q0 * q0 + q3 * q3 - 0.5)

        q0_ = -0.5 * q1 * p - 0.5 * q2 * q - 0.5 * q3 * r
        q1_ = 0.5 * q0 * p - 0.5 * q3 * q + 0.5 * q2 * r
        q2_ = 0.5 * q3 * p + 0.5 * q0 * q - 0.5 * q1 * r
        q3_ = 0.5 * q1 * q - 0.5 * q2 * p + 0.5 * q0 * r

        p_ = (
            -40.00063258437631 * pI - 2.8283979829540325 * p
        ) - 1.133407423682400 * q * r
        q_ = (
            -39.99980452524146 * qI - 2.8283752541008109 * q
        ) + 1.132078179613602 * p * r
        r_ = (
            -39.99978909742505 * rI - 2.8284134223281210 * r
        ) - 0.004695219977601 * p * q

        pI_ = p - pr
        qI_ = q - qr
        rI_ = r - rr
        hI_ = h - hr

        system_dynamics = [
            pn_,
            pe_,
            h_,
            u_,
            v_,
            w_,
            p_,
            q_,
            r_,
            q0_,
            q1_,
            q2_,
            q3_,
            pI_,
            qI_,
            rI_,
            hI_,
        ]  # has to be in the same order as the input variables
        # ===================================================

        return np.array(system_dynamics)  # return as numpy array


# CTRNN Fixed-Point Attractor example
# from https://easychair.org/publications/open/K6SZ
class CTRNN_fpa:
    def __init__(self, radius):
        # ============ adapt initial values ===========
        self.cx = (0.21535, -0.58587, 0.8, 0.52323, 0.5)  # initial values
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 0.01  # initial radius
        # ===================================================

        self.cx = np.array(self.cx, dtype=float)
        self.dim = self.cx.size  # dimension of the system

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim, dtype=object)

        # ============ adapt input and system dynamics ===========
        x1, x2, x3, x4, x5 = x

        x1p = (
            5419046323626097
            * (exp(-2 * x3) - 1)
            / (4503599627370496 * (exp(-2 * x3) + 1))
            - x1 / 1000000
            + 3601 * (exp(-2 * x4) - 1) / (50000 * (exp(-2 * x4) + 1))
            + 18727 * (exp(-2 * x5) - 1) / (20000 * (exp(-2 * x5) + 1))
        )
        x2p = (
            30003 * (exp(-2 * x4) - 1) / (20000 * (exp(-2 * x4) + 1))
            - 11881 * (exp(-2 * x3) - 1) / (10000 * (exp(-2 * x3) + 1))
            - x2 / 1000000
            - 93519 * (exp(-2 * x5) - 1) / (100000 * (exp(-2 * x5) + 1))
        )
        x3p = (
            7144123746377831
            * (exp(-2 * x3) - 1)
            / (4503599627370496 * (exp(-2 * x3) + 1))
            - x3 / 1000000
            - 5048886837752751
            * (exp(-2 * x4) - 1)
            / (72057594037927936 * (exp(-2 * x4) + 1))
            + 5564385670244745
            * (exp(-2 * x5) - 1)
            / (4503599627370496 * (exp(-2 * x5) + 1))
        )
        x4p = (
            1348796766312415
            * (exp(-2 * x4) - 1)
            / (4503599627370496 * (exp(-2 * x4) + 1))
            - 3086507593514335
            * (exp(-2 * x3) - 1)
            / (36028797018963968 * (exp(-2 * x3) + 1))
            - x4 / 1000000
            - 2476184452153819
            * (exp(-2 * x5) - 1)
            / (36028797018963968 * (exp(-2 * x5) + 1))
        )
        x5p = (
            1523758031023695
            * (exp(-2 * x4) - 1)
            / (18014398509481984 * (exp(-2 * x4) + 1))
            - 8060407538855891
            * (exp(-2 * x3) - 1)
            / (4503599627370496 * (exp(-2 * x3) + 1))
            - x5 / 1000000
            - 3139112893264555
            * (exp(-2 * x5) - 1)
            / (2251799813685248 * (exp(-2 * x5) + 1))
        )

        system_dynamics = [x1p, x2p, x3p, x4p, x5p]
        # ===================================================

        return np.array(system_dynamics)  # return as numpy array


# 12-dimensional cartpole with CT-RNN neural network controller
class CartpoleCTRNN:
    def __init__(self, radius):
        # ============ adapt initial values ===========
        self.cx = (0, 0, 0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 1e-4
        # ===================================================
        self.cx = np.array(self.cx, dtype=float)
        self.dim = self.cx.size  # dimension of the system

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim, dtype=object)

        # ============ adapt input and system dynamics ===========
        x_00, x_01, x_02, x_03, h_00, h_01, h_02, h_03, h_04, h_05, h_06, h_07 = x

        h_00_prime = -h_00 + tanh(
            h_00 * -1.49394
            + h_01 * -0.61947
            + h_02 * 0.37393
            + h_03 * -0.63451
            + h_04 * -1.08420
            + h_05 * 2.57981
            + h_06 * -1.53850
            + h_07 * -1.64354
            + x_00 * -0.07445
            + x_01 * 0.08736
            + x_02 * 0.47684
            + x_03 * 0.56397
            + 1.81867
        )
        h_01_prime = -h_01 + tanh(
            h_00 * 1.38295
            + h_01 * -1.45811
            + h_02 * 1.01473
            + h_03 * -0.04578
            + h_04 * -0.13416
            + h_05 * -0.21970
            + h_06 * 0.41791
            + h_07 * -0.10833
            + x_00 * -1.70409
            + x_01 * 0.51560
            + x_02 * -0.71273
            + x_03 * -0.61720
            + 0.86009
        )
        h_02_prime = -h_02 + tanh(
            h_00 * 0.57055
            + h_01 * 0.05941
            + h_02 * -0.16993
            + h_03 * -0.69688
            + h_04 * -0.30939
            + h_05 * -1.31558
            + h_06 * 0.03316
            + h_07 * 2.35873
            + x_00 * 0.21849
            + x_01 * 1.42990
            + x_02 * 1.60666
            + x_03 * -0.66847
            + -1.56112
        )
        h_03_prime = -h_03 + tanh(
            h_00 * -1.53087
            + h_01 * -0.42779
            + h_02 * -0.02195
            + h_03 * 0.18007
            + h_04 * 1.54262
            + h_05 * -0.19275
            + h_06 * -0.64598
            + h_07 * -0.85840
            + x_00 * 0.10095
            + x_01 * 0.55012
            + x_02 * 1.51416
            + x_03 * 1.95952
            + 0.00020
        )
        h_04_prime = -h_04 + tanh(
            h_00 * -0.54493
            + h_01 * 0.41139
            + h_02 * -0.18790
            + h_03 * -0.14312
            + h_04 * -0.10144
            + h_05 * 2.38792
            + h_06 * -0.44134
            + h_07 * -2.06462
            + x_00 * -0.01260
            + x_01 * 0.84681
            + x_02 * 0.51889
            + x_03 * 0.35089
            + -0.06552
        )
        h_05_prime = -h_05 + tanh(
            h_00 * 1.08793
            + h_01 * -0.24559
            + h_02 * -0.97041
            + h_03 * -0.16874
            + h_04 * -0.27133
            + h_05 * 0.73513
            + h_06 * 0.73096
            + h_07 * -2.26244
            + x_00 * 0.80312
            + x_01 * -1.35635
            + x_02 * -0.95283
            + x_03 * -0.68460
            + 0.11646
        )
        h_06_prime = -h_06 + tanh(
            h_00 * 0.08745
            + h_01 * 0.10104
            + h_02 * 0.82251
            + h_03 * -0.59823
            + h_04 * 1.16949
            + h_05 * 1.73511
            + h_06 * -0.89387
            + h_07 * 0.77149
            + x_00 * 0.22731
            + x_01 * -0.73872
            + x_02 * 2.83509
            + x_03 * -2.49535
            + 0.70343
        )
        h_07_prime = -h_07 + tanh(
            h_00 * -0.89045
            + h_01 * -1.30985
            + h_02 * 0.48740
            + h_03 * -0.45750
            + h_04 * -0.70727
            + h_05 * -0.43216
            + h_06 * -0.29159
            + h_07 * 3.70331
            + x_00 * -0.13294
            + x_01 * -0.23242
            + x_02 * 0.33244
            + x_03 * -1.45838
            + -0.02392
        )

        action_00 = tanh(
            h_00 * -0.32011
            + h_01 * 0.32958
            + h_02 * 0.07718
            + h_03 * 2.39431
            + h_04 * -1.37732
            + h_05 * 0.92722
            + h_06 * -1.07137
            + h_07 * -1.47286
        )

        gravity = 9.8
        masscart = 1
        masspole = 0.1
        total_mass = masscart + masspole
        length = 0.5
        polemass_length = masspole * length
        force_mag = 10

        force = force_mag * action_00
        costheta = cos(x_02)
        sintheta = sin(x_02)
        x_dot = x_01
        theta_dot = x_03

        temp = (force + polemass_length * theta_dot * theta_dot * sintheta) / total_mass
        thetaacc = (gravity * sintheta - costheta * temp) / (
            length * (4.0 / 3.0 - masspole * costheta * costheta / total_mass)
        )
        xacc = temp - polemass_length * thetaacc * costheta / total_mass

        x_00_prime = x_dot
        x_01_prime = xacc
        x_02_prime = theta_dot
        x_03_prime = thetaacc

        system_dynamics = [
            x_00_prime,
            x_01_prime,
            x_02_prime,
            x_03_prime,
            h_00_prime,
            h_01_prime,
            h_02_prime,
            h_03_prime,
            h_04_prime,
            h_05_prime,
            h_06_prime,
            h_07_prime,
        ]  # has to be in the same order as the input variables
        # ===================================================

        return np.array(system_dynamics)  # return as numpy array


# 12-dimensional cartpole with LTC neural network controller
class CartpoleLTC:
    def __init__(self, radius):
        # ============ adapt initial values ===========
        self.cx = (0, 0, 0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 1e-4
        # ===================================================
        self.cx = np.array(self.cx, dtype=float)
        self.dim = self.cx.size  # dimension of the system

    def true_sigmoid(self, x):
        return 0.5 * (tanh(x * 0.5) + 1)

    def sigmoid(self, v_pre, mu, sigma):
        mues = v_pre - mu
        x = sigma * mues
        return self.true_sigmoid(x)

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim, dtype=object)

        # ============ adapt input and system dynamics ===========
        x_00, x_01, x_02, x_03, h_00, h_01, h_02, h_03, h_04, h_05, h_06, h_07 = x

        u_00 = x_00 / 0.1
        u_01 = x_01 / 0.2
        u_02 = x_02 / 0.1
        u_03 = x_03 / 0.1

        swa_00 = (
            (-3.337809 - h_00) * 0.000100 * self.sigmoid(u_00, -1.822089, 4.920430)
            + (-3.793303 - h_00) * 0.217915 * self.sigmoid(u_01, 0.589409, 4.212314)
            + (1.074355 - h_00) * 0.000100 * self.sigmoid(u_02, -1.301882, 1.045275)
            + (-1.837527 - h_00) * 0.784721 * self.sigmoid(u_03, 0.061974, 9.342714)
        )
        swa_01 = (
            (1.922444 - h_01) * 0.048898 * self.sigmoid(u_00, 0.069332, 4.315088)
            + (-2.419416 - h_01) * 0.627589 * self.sigmoid(u_01, 2.746536, 3.680750)
            + (3.581162 - h_01) * 0.165919 * self.sigmoid(u_02, -1.045852, 3.066419)
            + (1.778340 - h_01) * 0.120764 * self.sigmoid(u_03, 1.837324, 3.448219)
        )
        swa_02 = (
            (2.205145 - h_02) * 0.304337 * self.sigmoid(u_00, 0.031854, 4.853640)
            + (-3.098629 - h_02) * 0.084116 * self.sigmoid(u_01, 0.441696, 3.296364)
            + (-1.226080 - h_02) * 0.177806 * self.sigmoid(u_02, 2.182991, 1.098078)
            + (-2.913366 - h_02) * 0.491636 * self.sigmoid(u_03, 2.224604, 1.728117)
        )
        swa_03 = (
            (0.960911 - h_03) * 0.364263 * self.sigmoid(u_00, -1.053052, 1.971431)
            + (-1.638120 - h_03) * 0.162840 * self.sigmoid(u_01, 7.642996, 2.541498)
            + (-1.138980 - h_03) * 0.746401 * self.sigmoid(u_02, 1.725237, 0.663583)
            + (1.642858 - h_03) * 0.299424 * self.sigmoid(u_03, 0.687038, -0.347229)
        )
        swa_04 = (
            (-0.752243 - h_04) * 0.649805 * self.sigmoid(u_00, 5.850876, 4.392767)
            + (0.135855 - h_04) * 0.018294 * self.sigmoid(u_01, 3.993045, 2.950581)
            + (-0.291421 - h_04) * 0.404132 * self.sigmoid(u_02, -0.929036, 3.571369)
            + (-3.123502 - h_04) * 0.521996 * self.sigmoid(u_03, -0.641402, 7.249665)
        )
        swa_05 = (
            (-0.472946 - h_05) * 0.359657 * self.sigmoid(u_00, 1.851295, 6.572009)
            + (-0.733541 - h_05) * 0.382768 * self.sigmoid(u_01, 2.120606, 3.201269)
            + (1.819771 - h_05) * 0.223484 * self.sigmoid(u_02, 2.934255, 4.131283)
            + (0.080200 - h_05) * 0.108734 * self.sigmoid(u_03, -2.828882, 1.771147)
        )
        swa_06 = (
            (2.596139 - h_06) * 0.233149 * self.sigmoid(u_00, 1.404586, 5.493755)
            + (-0.948662 - h_06) * 0.297568 * self.sigmoid(u_01, -5.682380, 3.978102)
            + (4.558042 - h_06) * 0.736431 * self.sigmoid(u_02, -0.067011, 3.915714)
            + (1.514289 - h_06) * 0.226179 * self.sigmoid(u_03, -1.101492, 5.705000)
        )
        swa_07 = (
            (0.907580 - h_07) * 0.077591 * self.sigmoid(u_00, 0.822831, 2.474382)
            + (-2.282039 - h_07) * 0.161999 * self.sigmoid(u_01, -0.241190, 6.074738)
            + (0.294892 - h_07) * 0.202664 * self.sigmoid(u_02, -0.289602, 8.727577)
            + (-2.329513 - h_07) * 0.242021 * self.sigmoid(u_03, -0.960229, 4.780694)
        )
        wa_00 = (
            (0.425183 - h_00) * 0.580120 * self.sigmoid(h_00, 1.888106, 5.698736)
            + (-0.026397 - h_00) * 0.422325 * self.sigmoid(h_01, 0.322996, 3.287380)
            + (0.232258 - h_00) * 0.234096 * self.sigmoid(h_02, 1.988594, 5.637971)
            + (0.325219 - h_00) * 0.822936 * self.sigmoid(h_03, 1.430856, 1.847167)
            + (-3.234136 - h_00) * 0.143390 * self.sigmoid(h_04, 1.693733, 1.529007)
            + (-0.446780 - h_00) * 0.000010 * self.sigmoid(h_05, 4.226530, 1.182250)
            + (-1.355023 - h_00) * 0.605004 * self.sigmoid(h_06, 1.083694, 2.011393)
            + (-3.184705 - h_00) * 0.338830 * self.sigmoid(h_07, 5.826199, 2.559444)
        )
        wa_01 = (
            (-1.902008 - h_01) * 0.151651 * self.sigmoid(h_00, 1.296603, 3.500569)
            + (-0.865412 - h_01) * 0.042350 * self.sigmoid(h_01, 1.769154, 10.179186)
            + (-3.485676 - h_01) * 0.259284 * self.sigmoid(h_02, 1.534717, 3.293725)
            + (-1.715340 - h_01) * 0.197946 * self.sigmoid(h_03, 0.264122, 5.787919)
            + (-3.960100 - h_01) * 0.346343 * self.sigmoid(h_04, 4.226974, 5.074825)
            + (-1.888838 - h_01) * 0.167991 * self.sigmoid(h_05, -0.192584, 4.463276)
            + (3.343740 - h_01) * 0.407932 * self.sigmoid(h_06, -2.019538, 2.096939)
            + (-3.828055 - h_01) * 0.392470 * self.sigmoid(h_07, -1.217940, 3.255837)
        )
        wa_02 = (
            (-0.149339 - h_02) * 0.274537 * self.sigmoid(h_00, -0.229752, 0.488301)
            + (0.411455 - h_02) * 0.026122 * self.sigmoid(h_01, 4.011376, 4.030908)
            + (3.896320 - h_02) * 0.157079 * self.sigmoid(h_02, 2.459823, 3.892246)
            + (0.530434 - h_02) * 0.059498 * self.sigmoid(h_03, 1.132183, 1.906741)
            + (-0.194996 - h_02) * 0.087192 * self.sigmoid(h_04, 0.000571, 0.949683)
            + (-0.407292 - h_02) * 0.194082 * self.sigmoid(h_05, -0.997799, 1.219973)
            + (-6.620962 - h_02) * 0.398726 * self.sigmoid(h_06, 0.469706, 5.065333)
            + (1.787329 - h_02) * 1.050713 * self.sigmoid(h_07, 4.342476, -0.298296)
        )
        wa_03 = (
            (-5.372879 - h_03) * 0.027963 * self.sigmoid(h_00, 2.814907, 7.856180)
            + (-0.646331 - h_03) * 0.149699 * self.sigmoid(h_01, 0.987782, 2.744861)
            + (0.116237 - h_03) * 0.135569 * self.sigmoid(h_02, 4.085198, 3.937947)
            + (-2.523393 - h_03) * 0.165372 * self.sigmoid(h_03, 0.703319, 5.054082)
            + (0.815988 - h_03) * 0.189418 * self.sigmoid(h_04, 3.352685, 7.961621)
            + (-1.706797 - h_03) * 0.502507 * self.sigmoid(h_05, -0.760482, 5.458648)
            + (-2.023418 - h_03) * 0.026645 * self.sigmoid(h_06, 4.557816, 2.245746)
            + (-0.889180 - h_03) * 0.468236 * self.sigmoid(h_07, -2.680705, 7.208436)
        )
        wa_04 = (
            (1.343493 - h_04) * 0.603259 * self.sigmoid(h_00, -0.901840, 3.872216)
            + (0.831619 - h_04) * 0.426593 * self.sigmoid(h_01, -3.089995, 5.196197)
            + (0.421409 - h_04) * 0.068839 * self.sigmoid(h_02, 3.991761, 0.784316)
            + (2.355011 - h_04) * 0.018658 * self.sigmoid(h_03, -0.024262, 7.865137)
            + (0.593414 - h_04) * 0.252146 * self.sigmoid(h_04, -1.463303, 3.132787)
            + (-0.637192 - h_04) * 0.118786 * self.sigmoid(h_05, -1.024096, 4.974855)
            + (-0.732044 - h_04) * 0.588030 * self.sigmoid(h_06, 1.033544, 6.298468)
            + (-0.194411 - h_04) * 0.679828 * self.sigmoid(h_07, 0.410277, 8.185476)
        )
        wa_05 = (
            (0.495678 - h_05) * 0.338934 * self.sigmoid(h_00, 1.896829, 3.323362)
            + (0.680690 - h_05) * 0.558593 * self.sigmoid(h_01, -0.032828, 3.472205)
            + (1.564447 - h_05) * 0.416106 * self.sigmoid(h_02, -0.448733, 4.371096)
            + (0.896225 - h_05) * 0.102958 * self.sigmoid(h_03, -1.402512, 7.281286)
            + (3.866403 - h_05) * 0.388992 * self.sigmoid(h_04, -1.392615, 1.017738)
            + (-1.880220 - h_05) * 0.075561 * self.sigmoid(h_05, 0.562915, 6.004653)
            + (-5.541122 - h_05) * 0.279411 * self.sigmoid(h_06, -0.679593, 2.001722)
            + (-0.065618 - h_05) * 0.475569 * self.sigmoid(h_07, -1.229657, 3.177203)
        )
        wa_06 = (
            (-1.465561 - h_06) * 0.122170 * self.sigmoid(h_00, 0.706757, 2.401071)
            + (-2.485563 - h_06) * 0.091596 * self.sigmoid(h_01, 0.335557, 3.252133)
            + (1.158366 - h_06) * 0.519633 * self.sigmoid(h_02, -1.169937, 3.588379)
            + (1.955986 - h_06) * 0.252199 * self.sigmoid(h_03, -1.974608, 2.870511)
            + (0.743124 - h_06) * 0.651011 * self.sigmoid(h_04, -0.210935, 7.228260)
            + (3.452268 - h_06) * 0.019733 * self.sigmoid(h_05, 3.515318, 3.879946)
            + (0.604325 - h_06) * 0.191627 * self.sigmoid(h_06, -3.070420, 4.221974)
            + (2.465013 - h_06) * 0.143574 * self.sigmoid(h_07, -0.961958, 3.002333)
        )
        wa_07 = (
            (-2.135510 - h_07) * 0.082957 * self.sigmoid(h_00, 0.857560, 4.560907)
            + (-0.719410 - h_07) * 0.121560 * self.sigmoid(h_01, 2.004791, 7.674825)
            + (-0.095160 - h_07) * 0.085265 * self.sigmoid(h_02, 0.915076, 3.384454)
            + (-6.402920 - h_07) * 0.017902 * self.sigmoid(h_03, 1.900790, 3.606339)
            + (0.911732 - h_07) * 0.018167 * self.sigmoid(h_04, 0.001252, 7.378105)
            + (-0.876065 - h_07) * 0.119567 * self.sigmoid(h_05, 2.130298, 7.560087)
            + (-1.520981 - h_07) * 0.012488 * self.sigmoid(h_06, -0.275924, 2.646271)
            + (0.734387 - h_07) * 0.088953 * self.sigmoid(h_07, -1.333799, 3.340692)
        )
        h_00_prime = 10000000.000000 * (0.663909 * (4.880308 - h_00) + swa_00 + wa_00)
        h_01_prime = 0.897653 * (2.666445 * (-0.506661 - h_01) + swa_01 + wa_01)
        h_02_prime = 0.513291 * (4.721010 * (-1.314806 - h_02) + swa_02 + wa_02)
        h_03_prime = 0.301197 * (0.073180 * (0.044130 - h_03) + swa_03 + wa_03)
        h_04_prime = 0.404649 * (1.368898 * (-1.209986 - h_04) + swa_04 + wa_04)
        h_05_prime = 1.131324 * (0.000010 * (-2.630356 - h_05) + swa_05 + wa_05)
        h_06_prime = 10000000.000000 * (0.539536 * (-0.529881 - h_06) + swa_06 + wa_06)
        h_07_prime = 0.525138 * (0.507019 * (0.567584 - h_07) + swa_07 + wa_07)

        action_00 = tanh(
            h_00 * -0.895448
            + h_01 * -0.795233
            + h_02 * 0.209421
            + h_03 * 0.143533
            + h_04 * -0.084270
            + h_05 * -0.342232
            + h_06 * 0.333434
            + h_07 * -0.032306
            + 0.757913
        )

        gravity = 9.8
        masscart = 1
        masspole = 0.1
        total_mass = masscart + masspole
        length = 0.5
        polemass_length = masspole * length
        force_mag = 10

        force = force_mag * action_00
        costheta = cos(x_02)
        sintheta = sin(x_02)
        x_dot = x_01
        theta_dot = x_03

        temp = (force + polemass_length * theta_dot * theta_dot * sintheta) / total_mass
        thetaacc = (gravity * sintheta - costheta * temp) / (
            length * (4.0 / 3.0 - masspole * costheta * costheta / total_mass)
        )
        xacc = temp - polemass_length * thetaacc * costheta / total_mass

        x_00_prime = x_dot
        x_01_prime = xacc
        x_02_prime = theta_dot
        x_03_prime = thetaacc

        system_dynamics = [
            x_00_prime,
            x_01_prime,
            x_02_prime,
            x_03_prime,
            h_00_prime,
            h_01_prime,
            h_02_prime,
            h_03_prime,
            h_04_prime,
            h_05_prime,
            h_06_prime,
            h_07_prime,
        ]  # has to be in the same order as the input variables
        # ===================================================

        return np.array(system_dynamics)  # return as numpy array


# 12-dimensional cartpole with LTC neural network controller (trained with RK integrator)
class CartpoleLTC_RK:
    def __init__(self, radius):
        # ============ adapt initial values ===========
        self.cx = (0, 0, 0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 1e-4
        # ===================================================
        self.cx = np.array(self.cx, dtype=float)
        self.dim = self.cx.size  # dimension of the system

    def true_sigmoid(self, x):
        return 0.5 * (tanh(x * 0.5) + 1)

    def sigmoid(self, v_pre, mu, sigma):
        mues = v_pre - mu
        x = sigma * mues
        return self.true_sigmoid(x)

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim, dtype=object)

        # ============ adapt input and system dynamics ===========
        x_00, x_01, x_02, x_03, h_00, h_01, h_02, h_03, h_04, h_05, h_06, h_07 = x

        u_00 = x_00 / 0.1
        u_01 = x_01 / 0.2
        u_02 = x_02 / 0.1
        u_03 = x_03 / 0.1

        swa_00 = (
            (-0.431604 - h_00)
            * 0.072614
            * 0.5
            * (tanh(0.5 * (u_00 - (0.487739)) * 5.650091) + 1)
            + (-1.100827 - h_00)
            * 0.581228
            * 0.5
            * (tanh(0.5 * (u_01 - (1.385196)) * 3.634596) + 1)
            + (0.174867 - h_00)
            * 0.111020
            * 0.5
            * (tanh(0.5 * (u_02 - (-0.993259)) * 5.331078) + 1)
            + (-0.510866 - h_00)
            * 0.189538
            * 0.5
            * (tanh(0.5 * (u_03 - (-1.174057)) * 5.951439) + 1)
        )
        swa_01 = (
            (-0.596286 - h_01)
            * 0.304088
            * 0.5
            * (tanh(0.5 * (u_00 - (-0.441833)) * 3.592494) + 1)
            + (-1.364699 - h_01)
            * 0.266900
            * 0.5
            * (tanh(0.5 * (u_01 - (-0.104207)) * 4.467682) + 1)
            + (0.723861 - h_01)
            * 0.282621
            * 0.5
            * (tanh(0.5 * (u_02 - (0.735169)) * 4.735680) + 1)
            + (-2.003241 - h_01)
            * 0.162848
            * 0.5
            * (tanh(0.5 * (u_03 - (1.824110)) * 4.624085) + 1)
        )
        swa_02 = (
            (-2.213958 - h_02)
            * 0.050379
            * 0.5
            * (tanh(0.5 * (u_00 - (1.028516)) * 4.883958) + 1)
            + (1.837608 - h_02)
            * 0.000100
            * 0.5
            * (tanh(0.5 * (u_01 - (0.476658)) * 3.920279) + 1)
            + (-0.869022 - h_02)
            * 0.000100
            * 0.5
            * (tanh(0.5 * (u_02 - (0.146460)) * 4.128851) + 1)
            + (0.970412 - h_02)
            * 0.142598
            * 0.5
            * (tanh(0.5 * (u_03 - (-0.057289)) * 4.129013) + 1)
        )
        swa_03 = (
            (3.316434 - h_03)
            * 0.000100
            * 0.5
            * (tanh(0.5 * (u_00 - (-0.281031)) * 4.200078) + 1)
            + (-0.286490 - h_03)
            * 0.136481
            * 0.5
            * (tanh(0.5 * (u_01 - (2.742276)) * 3.584777) + 1)
            + (-0.305795 - h_03)
            * 0.186325
            * 0.5
            * (tanh(0.5 * (u_02 - (1.382006)) * 4.178387) + 1)
            + (2.739162 - h_03)
            * 0.059435
            * 0.5
            * (tanh(0.5 * (u_03 - (-0.377228)) * 3.719664) + 1)
        )
        swa_04 = (
            (-1.856016 - h_04)
            * 0.031760
            * 0.5
            * (tanh(0.5 * (u_00 - (-0.689140)) * 3.121773) + 1)
            + (-1.060924 - h_04)
            * 0.233160
            * 0.5
            * (tanh(0.5 * (u_01 - (0.130687)) * 2.643270) + 1)
            + (0.459155 - h_04)
            * 0.143452
            * 0.5
            * (tanh(0.5 * (u_02 - (0.574488)) * 1.886950) + 1)
            + (-0.358857 - h_04)
            * 0.100690
            * 0.5
            * (tanh(0.5 * (u_03 - (1.888403)) * 4.933154) + 1)
        )
        swa_05 = (
            (1.980404 - h_05)
            * 0.323992
            * 0.5
            * (tanh(0.5 * (u_00 - (1.983044)) * 4.761758) + 1)
            + (-1.479298 - h_05)
            * 0.058762
            * 0.5
            * (tanh(0.5 * (u_01 - (0.921888)) * 6.088301) + 1)
            + (1.952959 - h_05)
            * 0.088192
            * 0.5
            * (tanh(0.5 * (u_02 - (0.347450)) * 4.452416) + 1)
            + (-1.889391 - h_05)
            * 0.035840
            * 0.5
            * (tanh(0.5 * (u_03 - (-1.054210)) * 5.084546) + 1)
        )
        swa_06 = (
            (-0.891564 - h_06)
            * 0.238208
            * 0.5
            * (tanh(0.5 * (u_00 - (0.668283)) * 4.712173) + 1)
            + (-0.434193 - h_06)
            * 0.102071
            * 0.5
            * (tanh(0.5 * (u_01 - (-0.555892)) * 3.985909) + 1)
            + (-1.002575 - h_06)
            * 0.211428
            * 0.5
            * (tanh(0.5 * (u_02 - (1.759687)) * 5.538335) + 1)
            + (0.961403 - h_06)
            * 0.000100
            * 0.5
            * (tanh(0.5 * (u_03 - (-0.823010)) * 5.778478) + 1)
        )
        swa_07 = (
            (-2.280897 - h_07)
            * 0.055498
            * 0.5
            * (tanh(0.5 * (u_00 - (0.378455)) * 3.488947) + 1)
            + (3.006480 - h_07)
            * 0.026770
            * 0.5
            * (tanh(0.5 * (u_01 - (-0.096777)) * 2.703730) + 1)
            + (0.428589 - h_07)
            * 0.131216
            * 0.5
            * (tanh(0.5 * (u_02 - (-0.058215)) * 4.788692) + 1)
            + (-0.848094 - h_07)
            * 0.156781
            * 0.5
            * (tanh(0.5 * (u_03 - (0.870154)) * 5.505928) + 1)
        )
        wa_00 = (
            (-2.061894 - h_00)
            * 0.255103
            * 0.5
            * (tanh(0.5 * (h_00 - (-0.599924)) * 4.079573) + 1)
            + (-2.919805 - h_00)
            * 0.111550
            * 0.5
            * (tanh(0.5 * (h_01 - (-0.003773)) * 3.577720) + 1)
            + (-1.169591 - h_00)
            * 0.076621
            * 0.5
            * (tanh(0.5 * (h_02 - (0.268934)) * 3.482080) + 1)
            + (2.610878 - h_00)
            * 0.358292
            * 0.5
            * (tanh(0.5 * (h_03 - (-0.801209)) * 1.404748) + 1)
            + (0.464382 - h_00)
            * 0.396781
            * 0.5
            * (tanh(0.5 * (h_04 - (-0.353826)) * 6.078532) + 1)
            + (-1.108612 - h_00)
            * 0.039682
            * 0.5
            * (tanh(0.5 * (h_05 - (-0.083197)) * 4.031213) + 1)
            + (-1.782495 - h_00)
            * 0.192295
            * 0.5
            * (tanh(0.5 * (h_06 - (-0.135476)) * 3.405039) + 1)
            + (-0.745238 - h_00)
            * 0.224384
            * 0.5
            * (tanh(0.5 * (h_07 - (0.838400)) * 3.514007) + 1)
        )
        wa_01 = (
            (0.323076 - h_01)
            * 0.132399
            * 0.5
            * (tanh(0.5 * (h_00 - (1.430406)) * 3.536494) + 1)
            + (-1.075827 - h_01)
            * 0.183413
            * 0.5
            * (tanh(0.5 * (h_01 - (0.444819)) * 3.535979) + 1)
            + (0.727482 - h_01)
            * 0.179833
            * 0.5
            * (tanh(0.5 * (h_02 - (1.784908)) * 4.561027) + 1)
            + (0.138851 - h_01)
            * 0.154461
            * 0.5
            * (tanh(0.5 * (h_03 - (1.996656)) * 5.773379) + 1)
            + (-0.711774 - h_01)
            * 0.143956
            * 0.5
            * (tanh(0.5 * (h_04 - (-1.484266)) * 2.894240) + 1)
            + (-0.365611 - h_01)
            * 0.268626
            * 0.5
            * (tanh(0.5 * (h_05 - (0.714121)) * 1.408464) + 1)
            + (1.591881 - h_01)
            * 0.189146
            * 0.5
            * (tanh(0.5 * (h_06 - (2.912474)) * 3.563775) + 1)
            + (2.669420 - h_01)
            * 0.000010
            * 0.5
            * (tanh(0.5 * (h_07 - (0.507458)) * 3.459168) + 1)
        )
        wa_02 = (
            (2.356203 - h_02)
            * 0.072446
            * 0.5
            * (tanh(0.5 * (h_00 - (0.818119)) * 4.080008) + 1)
            + (0.534163 - h_02)
            * 0.645405
            * 0.5
            * (tanh(0.5 * (h_01 - (0.951815)) * 4.642598) + 1)
            + (0.577408 - h_02)
            * 0.038336
            * 0.5
            * (tanh(0.5 * (h_02 - (-0.348300)) * 2.439894) + 1)
            + (-1.214589 - h_02)
            * 0.310216
            * 0.5
            * (tanh(0.5 * (h_03 - (0.429035)) * 5.454509) + 1)
            + (-1.883640 - h_02)
            * 0.026713
            * 0.5
            * (tanh(0.5 * (h_04 - (1.463999)) * 3.258402) + 1)
            + (-1.458897 - h_02)
            * 0.054829
            * 0.5
            * (tanh(0.5 * (h_05 - (2.008786)) * 1.817449) + 1)
            + (0.503231 - h_02)
            * 0.215241
            * 0.5
            * (tanh(0.5 * (h_06 - (-0.121074)) * 4.956780) + 1)
            + (-0.185760 - h_02)
            * 0.202772
            * 0.5
            * (tanh(0.5 * (h_07 - (-0.237805)) * 3.949310) + 1)
        )
        wa_03 = (
            (2.190941 - h_03)
            * 0.264046
            * 0.5
            * (tanh(0.5 * (h_00 - (2.558471)) * 6.171836) + 1)
            + (0.049259 - h_03)
            * 0.016450
            * 0.5
            * (tanh(0.5 * (h_01 - (-0.471219)) * 3.276491) + 1)
            + (1.897136 - h_03)
            * 0.093556
            * 0.5
            * (tanh(0.5 * (h_02 - (-0.834410)) * 2.282021) + 1)
            + (-1.203204 - h_03)
            * 0.016330
            * 0.5
            * (tanh(0.5 * (h_03 - (-0.229405)) * 2.725281) + 1)
            + (2.328495 - h_03)
            * 0.032663
            * 0.5
            * (tanh(0.5 * (h_04 - (2.172752)) * 1.766576) + 1)
            + (1.529989 - h_03)
            * 0.128894
            * 0.5
            * (tanh(0.5 * (h_05 - (-0.905643)) * 4.200280) + 1)
            + (-0.033709 - h_03)
            * 0.168944
            * 0.5
            * (tanh(0.5 * (h_06 - (0.228327)) * 1.674565) + 1)
            + (-0.903095 - h_03)
            * 0.101751
            * 0.5
            * (tanh(0.5 * (h_07 - (0.474075)) * 3.122658) + 1)
        )
        wa_04 = (
            (1.356759 - h_04)
            * 0.207902
            * 0.5
            * (tanh(0.5 * (h_00 - (-0.832245)) * 5.076437) + 1)
            + (-1.713048 - h_04)
            * 0.253930
            * 0.5
            * (tanh(0.5 * (h_01 - (-0.361454)) * 6.056005) + 1)
            + (1.365151 - h_04)
            * 0.025809
            * 0.5
            * (tanh(0.5 * (h_02 - (2.691158)) * 2.098762) + 1)
            + (-1.825929 - h_04)
            * 0.039653
            * 0.5
            * (tanh(0.5 * (h_03 - (-1.585948)) * 3.676085) + 1)
            + (-0.728560 - h_04)
            * 0.012267
            * 0.5
            * (tanh(0.5 * (h_04 - (-1.516692)) * 5.575822) + 1)
            + (-1.497486 - h_04)
            * 0.037167
            * 0.5
            * (tanh(0.5 * (h_05 - (0.021516)) * 5.903266) + 1)
            + (-2.542053 - h_04)
            * 0.181942
            * 0.5
            * (tanh(0.5 * (h_06 - (0.198382)) * 3.954360) + 1)
            + (3.068608 - h_04)
            * 0.095478
            * 0.5
            * (tanh(0.5 * (h_07 - (0.247832)) * 4.194288) + 1)
        )
        wa_05 = (
            (0.924399 - h_05)
            * 0.000010
            * 0.5
            * (tanh(0.5 * (h_00 - (1.051798)) * 4.877949) + 1)
            + (0.179591 - h_05)
            * 0.082066
            * 0.5
            * (tanh(0.5 * (h_01 - (0.472982)) * 5.943190) + 1)
            + (0.422308 - h_05)
            * 0.361843
            * 0.5
            * (tanh(0.5 * (h_02 - (0.984128)) * 4.885800) + 1)
            + (-1.484202 - h_05)
            * 0.660023
            * 0.5
            * (tanh(0.5 * (h_03 - (2.561464)) * 3.960099) + 1)
            + (3.444595 - h_05)
            * 0.269937
            * 0.5
            * (tanh(0.5 * (h_04 - (2.765615)) * 3.529421) + 1)
            + (-0.460907 - h_05)
            * 0.018725
            * 0.5
            * (tanh(0.5 * (h_05 - (0.951335)) * 7.133490) + 1)
            + (-1.438883 - h_05)
            * 0.135012
            * 0.5
            * (tanh(0.5 * (h_06 - (2.470155)) * 5.090559) + 1)
            + (-0.022624 - h_05)
            * 0.030127
            * 0.5
            * (tanh(0.5 * (h_07 - (-0.370460)) * 4.233265) + 1)
        )
        wa_06 = (
            (-2.204599 - h_06)
            * 0.144977
            * 0.5
            * (tanh(0.5 * (h_00 - (2.253059)) * 3.878229) + 1)
            + (0.148456 - h_06)
            * 0.270476
            * 0.5
            * (tanh(0.5 * (h_01 - (2.298548)) * 5.617167) + 1)
            + (0.544189 - h_06)
            * 0.187231
            * 0.5
            * (tanh(0.5 * (h_02 - (1.738353)) * 3.412084) + 1)
            + (0.957022 - h_06)
            * 0.091821
            * 0.5
            * (tanh(0.5 * (h_03 - (2.433656)) * 1.501055) + 1)
            + (0.305039 - h_06)
            * 0.000010
            * 0.5
            * (tanh(0.5 * (h_04 - (2.545896)) * 2.064504) + 1)
            + (2.075277 - h_06)
            * 0.172349
            * 0.5
            * (tanh(0.5 * (h_05 - (0.215808)) * 4.761319) + 1)
            + (0.651300 - h_06)
            * 0.392033
            * 0.5
            * (tanh(0.5 * (h_06 - (-2.850508)) * 5.109767) + 1)
            + (0.305318 - h_06)
            * 0.056583
            * 0.5
            * (tanh(0.5 * (h_07 - (-0.545839)) * 5.304969) + 1)
        )
        wa_07 = (
            (-2.419794 - h_07)
            * 0.164188
            * 0.5
            * (tanh(0.5 * (h_00 - (0.646132)) * 5.431365) + 1)
            + (-1.410748 - h_07)
            * 0.128581
            * 0.5
            * (tanh(0.5 * (h_01 - (1.941283)) * 2.716641) + 1)
            + (0.566556 - h_07)
            * 0.189314
            * 0.5
            * (tanh(0.5 * (h_02 - (2.011078)) * 6.978758) + 1)
            + (-0.999432 - h_07)
            * 0.038773
            * 0.5
            * (tanh(0.5 * (h_03 - (2.175668)) * 2.841975) + 1)
            + (-1.244396 - h_07)
            * 0.024732
            * 0.5
            * (tanh(0.5 * (h_04 - (-1.121634)) * 4.267510) + 1)
            + (0.211258 - h_07)
            * 0.069162
            * 0.5
            * (tanh(0.5 * (h_05 - (0.894198)) * 2.127854) + 1)
            + (2.880166 - h_07)
            * 0.023222
            * 0.5
            * (tanh(0.5 * (h_06 - (0.345695)) * 3.201268) + 1)
            + (-2.649584 - h_07)
            * 0.348009
            * 0.5
            * (tanh(0.5 * (h_07 - (0.311026)) * 3.550425) + 1)
        )
        h_00_prime = 0.472549 * (0.026911 * (-0.601392 - h_00) + swa_00 + wa_00)
        h_01_prime = 0.422074 * (0.001000 * (1.909540 - h_01) + swa_01 + wa_01)
        h_02_prime = 0.316776 * (1.660939 * (-0.450735 - h_02) + swa_02 + wa_02)
        h_03_prime = 0.484111 * (0.144001 * (1.478572 - h_03) + swa_03 + wa_03)
        h_04_prime = 0.414692 * (1.029787 * (-1.105074 - h_04) + swa_04 + wa_04)
        h_05_prime = 0.358391 * (0.218205 * (-1.095135 - h_05) + swa_05 + wa_05)
        h_06_prime = 0.880162 * (0.268332 * (0.043464 - h_06) + swa_06 + wa_06)
        h_07_prime = 1.671041 * (0.153414 * (0.695311 - h_07) + swa_07 + wa_07)

        action_00 = tanh(
            h_00 * 0.135081
            + h_01 * 0.189091
            + h_02 * -0.350403
            + h_03 * 0.388546
            + h_04 * 0.397087
            + h_05 * 0.186882
            + h_06 * -0.047590
            + h_07 * 0.197737
            + -0.352130
        )

        gravity = 9.8
        masscart = 1
        masspole = 0.1
        total_mass = masscart + masspole
        length = 0.5
        polemass_length = masspole * length
        force_mag = 10

        force = force_mag * action_00
        costheta = cos(x_02)
        sintheta = sin(x_02)
        x_dot = x_01
        theta_dot = x_03

        temp = (force + polemass_length * theta_dot * theta_dot * sintheta) / total_mass
        thetaacc = (gravity * sintheta - costheta * temp) / (
            length * (4.0 / 3.0 - masspole * costheta * costheta / total_mass)
        )
        xacc = temp - polemass_length * thetaacc * costheta / total_mass

        x_00_prime = x_dot
        x_01_prime = xacc
        x_02_prime = theta_dot
        x_03_prime = thetaacc

        system_dynamics = [
            x_00_prime,
            x_01_prime,
            x_02_prime,
            x_03_prime,
            h_00_prime,
            h_01_prime,
            h_02_prime,
            h_03_prime,
            h_04_prime,
            h_05_prime,
            h_06_prime,
            h_07_prime,
        ]  # has to be in the same order as the input variables
        # ===================================================

        return np.array(system_dynamics)  # return as numpy array


class TestNODE:
    def __init__(self, radius):
        # ============ adapt initial values ===========
        self.cx = (0, 0)
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 1e-4
        # ===================================================
        self.cx = np.array(self.cx, dtype=float)
        self.dim = self.cx.size  # dimension of the system

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim, dtype=object)

        # ============ adapt input and system dynamics ===========
        x0, x1 = x
        x0 = 0.5889 * tanh(0.4256 * x0 + 0.5061 * x1 + 0.1773) - 0.1000 * x0
        x1 = 0.3857 * tanh(-0.5563 * x0 + -0.1262 * x1 + -0.2136) - 0.1000 * x1
        system_dynamics = [x0, x1]
        # ===================================================
        return np.array(system_dynamics)  # return as numpy array


class LDSwithCTRNN:
    def __init__(self, radius):
        # ============ adapt initial values ===========
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 0.5
        # ===================================================
        self.cx = np.zeros(10)
        self.dim = self.cx.size  # dimension of the system
        arr = np.load("rl/lds_ctrnn.npz")
        self.params = {k: arr[k] for k in arr.files}

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim, dtype=object)

        hidden = np.tanh(np.dot(x, self.params["w1"]) + self.params["b1"])
        dhdt = np.dot(hidden, self.params["w2"]) + self.params["b2"]

        action = np.tanh(np.dot(hidden, self.params["wa"]) + self.params["ba"])
        x, y = x[-2:]

        dxdt = y
        dydt = 0.2 + 0.4 * action

        dxdt = np.array([dxdt]).reshape((1,))
        dydt = np.array([dydt]).reshape((1,))
        dfdt = np.concatenate([dhdt, dxdt, dydt], axis=0)
        return dfdt


class PendulumwithCTRNN:
    def __init__(self, radius):
        # ============ adapt initial values ===========
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 0.5
        # ===================================================
        self.cx = np.zeros(10)
        self.dim = self.cx.size  # dimension of the system
        arr = np.load("rl/pendulum_ctrnn.npz")
        self.params = {k: arr[k] for k in arr.files}

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim, dtype=object)

        hidden = np.tanh(np.dot(x, self.params["w1"]) + self.params["b1"])
        dhdt = np.dot(hidden, self.params["w2"]) + self.params["b2"]

        action = np.tanh(np.dot(hidden, self.params["wa"]) + self.params["ba"])
        th, thdot = x[-2:]

        max_speed = 8
        g = 9.81
        m = 1.0
        l = 1.0

        newthdot = -3 * g / (2 * l) * np.sin(th + np.pi) + 3.0 / (m * l ** 2) * action
        newthdot = max_speed * np.tanh(newthdot / max_speed)
        newth = newthdot

        dxdt = np.array([newth]).reshape((1,))
        dydt = np.array([newthdot]).reshape((1,))
        dfdt = np.concatenate([dhdt, dxdt, dydt], axis=0)
        return dfdt


class CTRNNosc:
    def __init__(self, radius):
        # ============ adapt initial values ===========
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 0.1
        # ===================================================
        self.cx = np.zeros(16)
        self.dim = self.cx.size  # dimension of the system
        arr = np.load("rl/ctrnn_osc.npz")
        self.params = {k: arr[k] for k in arr.files}

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim, dtype=object)

        hidden = np.tanh(np.dot(x, self.params["w1"]) + self.params["b1"])
        dhdt = np.dot(hidden, self.params["w2"]) + self.params["b2"]
        return dhdt

class spiralL:
    def __init__(self,radius):
        self.cx = (2.0, 0.0)
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 0.01
        self.cx = np.array(self.cx, dtype=float)
        self.dim = self.cx.size

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim, dtype=object)

        w1 = np.array([[-0.32294768,0.59955627],
        [0.47014388,-0.39748120],
        [-0.56326932,0.33752987],
        [0.45147443,0.31528524],
        [0.41403031,-0.47271276],
        [-0.12952870,-0.62095606],
        [-0.41343114,-0.45678866],
        [-0.33266136,0.29245856],
        [0.50114638,0.39612201],
        [0.47665390,0.55137879]])
        b1 = np.array([0.0038923009,0.0037905588,0.0017197595,-0.0033185149,0.0024190384,-0.0013056855,0.0011365928,-0.00042518601,-0.0025141449,0.0010660964])
        w2 = np.array([[-0.50525320,0.34800902,-0.34015974,-0.40054744,0.39193857,0.59363592,0.56743664,-0.33811751,-0.36945280,-0.46805024],
        [-0.41715327,0.56257814,-0.56921810,0.60423535,0.53992182,-0.14412111,-0.45906776,-0.35295558,0.49238238,0.43526673]])
        b2 = np.array([-0.0013696412,0.00060380378]);

        hidden = np.dot(x,np.transpose(w1)) + b1
        dhdt = np.dot(hidden, np.transpose(w2)) + b2
        return dhdt



class spiralNL:
    def __init__(self,radius):
        self.cx = (2.0, 0.0)
        if radius is not None:
            self.rad = radius
        else:
            self.rad = 0.01
        self.cx = np.array(self.cx, dtype=float)
        self.dim = self.cx.size

    def fdyn(self, t=0, x=None):
        if x is None:
            x = np.zeros(self.dim,dtype=object)
        w1 = np.array([[0.2911133 ,  0.12008807],
        [-0.24582624,  0.23181419],
        [-0.25797904,  0.21687193],
        [-0.19282854, -0.2602416],
        [0.26780415, -0.20697702],
        [0.23462369,  0.2294843],
        [-0.2583547 ,  0.21444395],
        [-0.04514714,  0.29514763],
        [-0.15318371, -0.275755],
        [0.24873598,  0.21018365]])
        b1 = np.array([0.0038677 , -0.00026365, -0.007168970,  0.02469357,  0.01338706,0.00856025, -0.00888401,  0.00516089, -0.00634514, -0.01914518])
        w2 = np.array([[-0.58693904, -0.814841  , -0.8175157 ,  0.97060364,  0.6908913 ,-0.92446184, -0.79249185, -1.1507587 ,  1.2072723 , -0.7983982],
        [1.1564877 , -0.8991244 , -1.0774536 , -0.6731967 ,  1.0154784 , 0.8984464 , -1.0766245 , -0.238209  , -0.5233613 ,  0.8886671]])
        b2 = np.array([-0.04129209, -0.01508532])

        hidden = np.tanh(np.dot(x,w1.T) + b1)
        dhdt = np.dot(hidden, w2.T) + b2
        return dhdt
