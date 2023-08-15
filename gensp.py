""" XMAP: Generate spdata file used by GTC and Orbit.

---------------
    Author: Xishuo Wei. (xishuow@uci.edu, weixishuo@gmail.com)
    
"""

import numpy as np
from fortran import spline_function_f
import logging
import pickle


def fortran_write(map_data, output):
    lsp = map_data["lsp"]
    lst = map_data["lst"]

    logging.info("construt bsp")
    bsp = np.zeros((9, lsp, lst), dtype="float32")
    bsp[0, :, :] = np.sqrt(map_data["BR_onptb"][:, :]**2 +
                           map_data["BZ_onptb"][:, :]**2 + map_data["BT_onptb"][:, :]**2)

    spdpsi = map_data["psimesh"][1]-map_data["psimesh"][0]
    spdtheta = map_data["tb_mesh"][1]-map_data["tb_mesh"][0]
    bsp = spline_function_f.spline_function.construct_spline2d(
        spdpsi, spdtheta, bsp, 1, 2, lsp, lst)

    logging.info("construt xsp")
    xsp = np.zeros((9, lsp, lst), dtype="float32")
    xsp[0, :, :] = map_data["Ronptb"][:, :]
    # xsp[0,:,:] = Ronptp[:,:]
    xsp = spline_function_f.spline_function.construct_spline2d(
        spdpsi, spdtheta, xsp, 1, 2, lsp, lst)

    logging.info("construt zsp")
    zsp = np.zeros((9, lsp, lst), dtype="float32")
    zsp[0, :, :] = map_data["Zonptb"][:, :] - map_data["Z_ctr"]
    zsp = spline_function_f.spline_function.construct_spline2d(
        spdpsi, spdtheta, zsp, 1, 2, lsp, lst)

    logging.info("construct delta")
    delsp = np.zeros((9, lsp, lst), dtype="float32")
    delsp[0, :, :] = map_data["delta"][:, :]
    delsp = spline_function_f.spline_function.construct_spline2d(
        spdpsi, spdtheta, delsp, 0, 2, lsp, lst)

    logging.info("construct fsp")
    fsp = np.zeros((9, lsp, lst), dtype="float32")
    fsp[0, :, :] = map_data["nu_neg"][:, :]
    fsp = spline_function_f.spline_function.construct_spline2d(
        spdpsi, spdtheta, fsp, 0, 2, lsp, lst)

    logging.info("construt jacsp")
    jacsp = np.zeros((9, lsp, lst), dtype="float32")
    for jt in range(0, lst):
        jacsp[0, :, jt] = (map_data["gpsi"] * map_data["qpsi"] +
                           map_data["cpsi"]) / bsp[0, :, jt]

    jacsp = spline_function_f.spline_function.construct_spline2d(
        spdpsi, spdtheta, jacsp, 1, 2, lsp, lst)

    logging.info("constuct rd2d")
    rd2d = np.zeros((9, lsp, lst), dtype="float32")
    for jt in range(0, lst):
        rd2d[0, :, jt] = map_data["cpsi"]
    rd2d = spline_function_f.spline_function.construct_spline2d(
        spdpsi, spdtheta, rd2d, 0, 2, lsp, lst)

    logging.info("construt qpsi")
    qpsi = np.zeros((3, lsp), dtype="float32")
    qpsi[0, :] = map_data["qpsi"]
    qpsi = spline_function_f.spline_function.construct_spline1d(
        spdpsi, qpsi, 0, lsp)

    logging.info("construct gpsi")
    gpsi = np.zeros((3, lsp), dtype="float32")
    gpsi[0, :] = map_data["gpsi"]
    gpsi = spline_function_f.spline_function.construct_spline1d(
        spdpsi, gpsi, 0, lsp)

    logging.info("construct ppsi")
    ppsi = np.zeros((3, lsp), dtype="float32")
    ppsi[0, :] = map_data["ppsi"]
    ppsi = spline_function_f.spline_function.construct_spline1d(
        spdpsi, ppsi, 0, lsp)

    logging.info("construct rpsi")
    rpsi = np.zeros((3, lsp), dtype="float32")
    rpsi[0, :] = (map_data["Ronptb"][:, 0]-map_data["R_ctr"]) / \
        map_data["R_ctr"]
    rpsi = spline_function_f.spline_function.construct_spline1d(
        spdpsi, rpsi, 1, lsp)

    logging.info("construct torpsi")
    torpsi = np.zeros((3, lsp), dtype="float32")
    for ip in range(1, lsp):
        torpsi[0, ip] = torpsi[0, ip-1] + qpsi[0, ip-1] * spdpsi + 0.5*qpsi[1, ip-1] * spdpsi**2\
            + 1/3*qpsi[2, ip-1] * spdpsi**3

    torpsi = spline_function_f.spline_function.construct_spline1d(
        spdpsi, torpsi, 0, lsp)

    # write toroidal psi for profile.dat
#    with open("torpsispickle", "wb") as fb:
#        pickle.dump(map_data["psimesh"], fb)
#        pickle.dump(torpsi, fb)

    # write to spdata

    from fortran import spdata_interface_f

    output_file = output
    logging.info(f"Writing spdata to {output_file}...")

    spdata_interface_f.spdataio.filename(output_file)
    spdata_interface_f.spdataio.allocate(lsp, lst)

    spdata_interface_f.spdataio.bsp = bsp
    spdata_interface_f.spdataio.xsp = xsp
    spdata_interface_f.spdataio.zsp = zsp
    spdata_interface_f.spdataio.gsp = jacsp
    spdata_interface_f.spdataio.qpsi = qpsi
    spdata_interface_f.spdataio.gpsi = gpsi
    spdata_interface_f.spdataio.rd = rd2d
    spdata_interface_f.spdataio.ppsi = ppsi
    spdata_interface_f.spdataio.rpsi = rpsi
    spdata_interface_f.spdataio.fsp = fsp
    spdata_interface_f.spdataio.delsp = delsp
    spdata_interface_f.spdataio.torpsi = torpsi
    spdata_interface_f.spdataio.psiw = map_data["psimesh"][-1]
    spdata_interface_f.spdataio.ped = map_data["psimesh"][-1]

    spdata_interface_f.spdataio.write()


def python_write(map_data, output):
    import fortranformat as ff
    lsp = map_data["lsp"]
    lst = map_data["lst"]
    logger = logging.getLogger('xmap')

    logger.info("construt bsp")
    bsp = np.zeros((9, lsp, lst), dtype="float32")
    bsp[0, :, :] = np.sqrt(map_data["BR_onptb"][:, :]**2 +
                           map_data["BZ_onptb"][:, :]**2 + map_data["BT_onptb"][:, :]**2)

    spdpsi = map_data["psimesh"][1]-map_data["psimesh"][0]
    spdtheta = map_data["tb_mesh"][1]-map_data["tb_mesh"][0]
    bsp = spline_function_f.spline_function.construct_spline2d(
        spdpsi, spdtheta, bsp, 1, 2, lsp, lst)

    logger.info("construt xsp")
    xsp = np.zeros((9, lsp, lst), dtype="float32")
    xsp[0, :, :] = map_data["Ronptb"][:, :]
    # xsp[0,:,:] = Ronptp[:,:]
    xsp = spline_function_f.spline_function.construct_spline2d(
        spdpsi, spdtheta, xsp, 1, 2, lsp, lst)

    logger.info("construt zsp")
    zsp = np.zeros((9, lsp, lst), dtype="float32")
    zsp[0, :, :] = map_data["Zonptb"][:, :] - map_data["Z_ctr"]
    # zsp[0,:,:] = Zonptp[:,:]
    zsp = spline_function_f.spline_function.construct_spline2d(
        spdpsi, spdtheta, zsp, 1, 2, lsp, lst)

    logger.info("construct delta")
    delsp = np.zeros((9, lsp, lst), dtype="float32")
    delsp[0, :, :] = map_data["delta"][:, :]
    delsp = spline_function_f.spline_function.construct_spline2d(
        spdpsi, spdtheta, delsp, 0, 2, lsp, lst)

    logger.info("construct fsp")
    fsp = np.zeros((9, lsp, lst), dtype="float32")
    fsp[0, :, :] = map_data["nu_neg"][:, :]
    fsp = spline_function_f.spline_function.construct_spline2d(
        spdpsi, spdtheta, fsp, 0, 2, lsp, lst)

    logger.info("construt jacsp")
    jacsp = np.zeros((9, lsp, lst), dtype="float32")
    for jt in range(0, lst):
        jacsp[0, :, jt] = (map_data["gpsi"] * map_data["qpsi"] +
                           map_data["cpsi"]) / bsp[0, :, jt]

    jacsp = spline_function_f.spline_function.construct_spline2d(
        spdpsi, spdtheta, jacsp, 1, 2, lsp, lst)

    logger.info("constuct rd2d")
    rd2d = np.zeros((9, lsp, lst), dtype="float32")
    for jt in range(0, lst):
        rd2d[0, :, jt] = map_data["cpsi"]
    rd2d = spline_function_f.spline_function.construct_spline2d(
        spdpsi, spdtheta, rd2d, 0, 2, lsp, lst)

    logger.info("construt qpsi")
    qpsi = np.zeros((3, lsp), dtype="float32")
    qpsi[0, :] = map_data["qpsi"]
    qpsi = spline_function_f.spline_function.construct_spline1d(
        spdpsi, qpsi, 0, lsp)

    logger.info("construct gpsi")
    gpsi = np.zeros((3, lsp), dtype="float32")
    gpsi[0, :] = map_data["gpsi"]
    gpsi = spline_function_f.spline_function.construct_spline1d(
        spdpsi, gpsi, 0, lsp)

    logger.info("construct ppsi")
    ppsi = np.zeros((3, lsp), dtype="float32")
    ppsi[0, :] = map_data["ppsi"]
    ppsi = spline_function_f.spline_function.construct_spline1d(
        spdpsi, ppsi, 0, lsp)

    logger.info("construct rpsi")
    rpsi = np.zeros((3, lsp), dtype="float32")
    rpsi[0, :] = (map_data["Ronptb"][:, 0]-map_data["R_ctr"]) / \
        map_data["R_ctr"]
    rpsi = spline_function_f.spline_function.construct_spline1d(
        spdpsi, rpsi, 1, lsp)

    logger.info("construct torpsi")
    torpsi = np.zeros((3, lsp), dtype="float32")
    for ip in range(1, lsp):
        torpsi[0, ip] = torpsi[0, ip-1] + qpsi[0, ip-1] * spdpsi + 0.5*qpsi[1, ip-1] * spdpsi**2\
            + 1/3*qpsi[2, ip-1] * spdpsi**3

    torpsi = spline_function_f.spline_function.construct_spline1d(
        spdpsi, torpsi, 0, lsp)
    logger.info(f"Writing spdata...")

    from datetime import date
    today = date.today()
    output.write(f"{today} Generated by XMAP. Author: weixishuo@gmail.com\n")
    int_line = ff.FortranRecordWriter("6I4")
    float_line = ff.FortranRecordWriter("1p4e18.10")
    output.write(int_line.write([lsp, lst-1, 0, 0]))
    output.write("\n")
    psiw = ped = map_data["psimesh"][-1]
    output.write(float_line.write([psiw, ped]))
    output.write("\n")
    for i in range(lsp):
        output.write(float_line.write(bsp[0, i, :]))  # B-field
        output.write("\n")
        output.write(float_line.write(bsp[1, i, :]))
        output.write("\n")
        output.write(float_line.write(bsp[2, i, :]))
        output.write("\n")
        output.write(float_line.write(bsp[3, i, :]))
        output.write("\n")
        output.write(float_line.write(bsp[4, i, :]))
        output.write("\n")
        output.write(float_line.write(bsp[5, i, :]))
        output.write("\n")
        output.write(float_line.write(bsp[6, i, :]))
        output.write("\n")
        output.write(float_line.write(bsp[7, i, :]))
        output.write("\n")
        output.write(float_line.write(bsp[8, i, :]))
        output.write("\n")

        output.write(float_line.write(xsp[0, i, :]))  # X-coordinate
        output.write("\n")
        output.write(float_line.write(xsp[1, i, :]))
        output.write("\n")
        output.write(float_line.write(xsp[2, i, :]))
        output.write("\n")
        output.write(float_line.write(xsp[3, i, :]))
        output.write("\n")
        output.write(float_line.write(xsp[4, i, :]))
        output.write("\n")
        output.write(float_line.write(xsp[5, i, :]))
        output.write("\n")
        output.write(float_line.write(xsp[6, i, :]))
        output.write("\n")
        output.write(float_line.write(xsp[7, i, :]))
        output.write("\n")
        output.write(float_line.write(xsp[8, i, :]))
        output.write("\n")

        output.write(float_line.write(zsp[0, i, :]))  # Z-coordinate
        output.write("\n")
        output.write(float_line.write(zsp[1, i, :]))
        output.write("\n")
        output.write(float_line.write(zsp[2, i, :]))
        output.write("\n")
        output.write(float_line.write(zsp[3, i, :]))
        output.write("\n")
        output.write(float_line.write(zsp[4, i, :]))
        output.write("\n")
        output.write(float_line.write(zsp[5, i, :]))
        output.write("\n")
        output.write(float_line.write(zsp[6, i, :]))
        output.write("\n")
        output.write(float_line.write(zsp[7, i, :]))
        output.write("\n")
        output.write(float_line.write(zsp[8, i, :]))
        output.write("\n")

        output.write(float_line.write(jacsp[0, i, :]))  # Gicobian
        output.write("\n")
        output.write(float_line.write(jacsp[1, i, :]))
        output.write("\n")
        output.write(float_line.write(jacsp[2, i, :]))
        output.write("\n")
        output.write(float_line.write(jacsp[3, i, :]))
        output.write("\n")
        output.write(float_line.write(jacsp[4, i, :]))
        output.write("\n")
        output.write(float_line.write(jacsp[5, i, :]))
        output.write("\n")
        output.write(float_line.write(jacsp[6, i, :]))
        output.write("\n")
        output.write(float_line.write(jacsp[7, i, :]))
        output.write("\n")
        output.write(float_line.write(jacsp[8, i, :]))
        output.write("\n")

        output.write(float_line.write(qpsi[:, i]))  # q
        output.write("\n")
        output.write(float_line.write(gpsi[:, i]))  # g-current
        output.write("\n")
        output.write(float_line.write(rd2d[0:3, i, 0]))  # I-current
        output.write("\n")
        output.write(float_line.write(ppsi[:, i]))  # pressure
        output.write("\n")
        output.write(float_line.write(rpsi[:, i]))  # minor radius
        output.write("\n")
        output.write(float_line.write(torpsi[:, i]))  # toroidal flux
        output.write("\n")

    output.write(int_line.write([0, 0]))
    output.write("\n")
    output.write(float_line.write([0.0, 0.0, 0.0]))
    output.write("\n")
    output.write(float_line.write([0.0, 0.0]))
    output.write("\n")

    for i in range(lsp):
        output.write(float_line.write(fsp[0, i, :]))  # Nu
        output.write("\n")
        output.write(float_line.write(fsp[1, i, :]))
        output.write("\n")
        output.write(float_line.write(fsp[2, i, :]))
        output.write("\n")
        output.write(float_line.write(fsp[3, i, :]))
        output.write("\n")
        output.write(float_line.write(fsp[4, i, :]))
        output.write("\n")
        output.write(float_line.write(fsp[5, i, :]))
        output.write("\n")
        output.write(float_line.write(fsp[6, i, :]))
        output.write("\n")
        output.write(float_line.write(fsp[7, i, :]))
        output.write("\n")
        output.write(float_line.write(fsp[8, i, :]))
        output.write("\n")

    for i in range(lsp):
        output.write(float_line.write(delsp[0, i, :]))  # delta-current
        output.write("\n")
        output.write(float_line.write(delsp[1, i, :]))
        output.write("\n")
        output.write(float_line.write(delsp[2, i, :]))
        output.write("\n")
        output.write(float_line.write(delsp[3, i, :]))
        output.write("\n")
        output.write(float_line.write(delsp[4, i, :]))
        output.write("\n")
        output.write(float_line.write(delsp[5, i, :]))
        output.write("\n")
        output.write(float_line.write(delsp[6, i, :]))
        output.write("\n")
        output.write(float_line.write(delsp[7, i, :]))
        output.write("\n")
        output.write(float_line.write(delsp[8, i, :]))
        output.write("\n")
