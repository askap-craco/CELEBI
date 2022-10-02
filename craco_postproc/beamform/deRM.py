def deRM(f, RM):
    freqs = np.arange(f.f0.value - 336 / 2, f.f0.value + 336 / 2) * u.MHz

    wlens = c.c / freqs
    lambda0 = c.c / f.f0

    psi = RM * (wlens ** 2 - lambda0 ** 2)
    psi_reshape = psi[:, np.newaxis].repeat(100000, axis=1)

    Q_deRM = f.q_ds * np.cos(2 * psi_reshape) + f.u_ds * np.sin(
        2 * psi_reshape
    )
    U_deRM = f.u_ds * np.cos(2 * psi_reshape) - f.q_ds * np.sin(
        2 * psi_reshape
    )

    return Q_deRM, U_deRM, psi_reshape
