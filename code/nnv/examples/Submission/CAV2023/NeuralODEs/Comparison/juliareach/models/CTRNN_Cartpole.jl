# ===================================
# Add description of benchmark
# Source
# ===================================

using ReachabilityAnalysis, Plots

@taylorize function CTRNN_Cartpole!(dx, x, p, t)
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

        dx[1] = x_00_prime
        dx[2] = x_01_prime
        dx[3] = x_02_prime
        dx[4] = x_03_prime
        dx[5] = h_00_prime
        dx[6] = h_01_prime
        dx[7] = h_02_prime
        dx[8] = h_03_prime
        dx[9] = h_04_prime
        dx[10] = h_05_prime
        dx[11] = h_06_prime
        dx[12] = h_07_prime
        # has to be in the same order as the input variables
        # ===================================================

	return dx
end

function CTRNN_Cartpole()
	# initial state
	X0 = (-0.0001 .. 0.0001) × (-0.0001 .. 0.0001) × (0.0009 .. 0.0011) × (-0.0001 .. 0.0001) × (-0.0001 .. 0.0001) × (-0.0001 .. 0.0001) × (-0.0001 .. 0.0001) × (-0.0001 .. 0.0001) × (-0.0001 .. 0.0001) × (-0.0001 .. 0.0001) × (-0.0001 .. 0.0001) × (-0.0001 .. 0.0001)
	return @ivp(x' = CTRNN_Cartpole!(x), dim:12, x(0) ∈ X0)
end
