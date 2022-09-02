""" XMAP: Generate spdata file used by GTC and Orbit.

---------------
    Author: Xishuo Wei. (xishuow@uci.edu, weixishuo@gmail.com)
    
"""


def check_plot(check_data, figs):
    import matplotlib
    if figs == "save":
        matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    psimesh = check_data["psimesh"]
    smooth_q = check_data["q1d"]
    mean_q2d = check_data["mean_q2d"]
    var_q2d = check_data["var_q2d"]

    fig = plt.figure(figsize=(8, 6))
    plt.plot(psimesh, smooth_q, linewidth=2, color="tab:blue")
    plt.errorbar(psimesh, mean_q2d, var_q2d, color="tab:orange",
                 ecolor="r", capsize=3, errorevery=5)
    plt.legend(
        [r"$q$", r"$B\cdot{\nabla\zeta}/B\cdot{\nabla\theta}$"], fontsize=20)
    plt.xlabel(r"$\psi_p$")
    if figs == "save":
        plt.savefig("q_check.png")

    cpsi = check_data["ri1d"]
    mean_ri2d = check_data["mean_ri2d"]
    var_ri2d = check_data["var_ri2d"]
    fig = plt.figure(figsize=(8, 6))
    plt.plot(psimesh, cpsi, linewidth=2, color="tab:blue")
    plt.errorbar(psimesh, mean_ri2d, var_ri2d, color="tab:orange",
                 ecolor="r", capsize=3, errorevery=5)
    plt.legend([r"$\bar{I}$", r"$B\cdot e_\theta$"], fontsize=20)
    plt.xlabel(r"$\psi_p$")
    if figs == "save":
        plt.savefig("I_check.png")

    gqi = check_data["gqi"]
    mean_Jacb2 = check_data["mean_Jacb2"]
    var_Jacb2 = check_data["var_Jacb2"]
    fig = plt.figure(figsize=(8, 6))
    plt.plot(psimesh, gqi, color="tab:blue")
    plt.errorbar(psimesh[1:], mean_Jacb2[1:], var_Jacb2,
                 color="tab:orange", ecolor="r", capsize=3, errorevery=5)
    plt.legend(["(gq+I)", "J*B^2"], fontsize=20)
    plt.xlabel(r"$\psi_p$")
    if figs == "save":
        plt.savefig("J_check.png")

    Ronptb = check_data["R2d"]
    Zonptb = check_data["Z2d"]
    rbdry = check_data["Rbdry"]
    zbdry = check_data["Zbdry"]
    xlim = check_data["xlim"]
    ylim = check_data["ylim"]
    npsi = check_data["npsi"]
    ntheta = check_data["ntheta"]

    fig, ax = plt.subplots(dpi=120)
    ax.set_aspect('equal', 'box')
    ax.plot(xlim, ylim, 'k-', linewidth=0.5)
    psi_end = npsi+1

    for jt in range(0, ntheta+1, 5):
        plt.plot(Ronptb[:psi_end, jt], Zonptb[:psi_end, jt], 'b-')
    for ip in range(0, psi_end, 10):
        plt.plot(Ronptb[ip, :], Zonptb[ip, :], 'k-')
    ax.plot(rbdry, zbdry, 'r--')
    if figs == "save":
        plt.savefig("Boozer.png")

    delta = check_data["del2d"]
    nu_onptb = check_data["nu2d"]

    fig, ax = plt.subplots(dpi=120)
    ax.set_aspect('equal', 'box')
    ax.pcolormesh(Ronptb, Zonptb, nu_onptb[:-1, :-1])
    ax.set_title(r"$\nu$")
    if figs == "save":
        plt.savefig("nu.png")

    fig, ax = plt.subplots(dpi=120)
    ax.set_aspect('equal', 'box')
    ax.pcolormesh(Ronptb, Zonptb, delta[:-1, :-1])
    ax.set_title(r"$\delta$")
    if figs == "save":
        plt.savefig("delta.png")

    if figs == "show":
        plt.show()
