""" XMAP: Generate spdata file used by GTC and Orbit.

---------------
    Author: Xishuo Wei. (xishuow@uci.edu, weixishuo@gmail.com)
    
"""

import argparse
from mapping import mapping_core
import logging
from gensp import fortran_write
from consistency_plot import check_plot
from geqdsk import read
from utils import FigType

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ini", help="A file containing parameters."
                        "Not necessary if parameters inputs in command lines.")
    parser.add_argument("--input", help="Path of input gfile")
    parser.add_argument("--output", help="Path of output spdata file")
    parser.add_argument("--lsp", help="Grid# in psi direction", type=int)
    parser.add_argument(
        "--lst", help="Grid# in theta direction. Must be an even number", type=int)
    parser.add_argument(
        "--psimax", help="ratio of psimax to LCFS, default=0.995", type=float, default=0.995)
    parser.add_argument(
        "--detail-print", help="Set to False if no debug info is wanted. default: False", default=False, type=bool)
    parser.add_argument(
        "--figs", type=FigType, choices=list(FigType), help="Set to 'show' for popup plots. Set to 'save' to save figures. default: show", default="show")
    args = parser.parse_args()

    if args.ini:
        import json
        import sys
        with open(args.ini, "r") as fj:
            paras = json.load(fj)
            filenames = paras["gfilename"]
            for ind, inputfile in enumerate(filenames):
                data, check_data = mapping_core(
                    inputfile, paras["npsi"], paras["ntheta"], paras["figs"])
                if paras["figs"]:
                    check_plot(check_data, paras["figs"])
                fortran_write(data)
            sys.exit()

    paras = args
    if paras.detail_print:
        logging.basicConfig(level=logging.INFO)
    if not paras.input:
        filename = input("Please input path of gfile:\n")
    else:
        filename = paras.input

    if not paras.lsp:
        lsp = input("Please input npsi (=lsp in GTC and spdata), default: 91\n")
        if not lsp:
            lsp = 91
        else:
            lsp = int(lsp)
    else:
        lsp = paras.lsp

    if not paras.lst:
        lst = input(
            "Please input ntheta (=lst in GTC, =lst+1 in spdata), default: 122\n")
        if not lst:
            lst = 122
        else:
            lst = int(lst)
    else:
        lst = paras.lst

    psimax = paras.psimax
    print(f"Ratio of last mapping flux surface to LCFS: {psimax}."
          "Use --psimax=[value] to specify other value.")

    if not paras.output:
        output = input(
            "Please enter the output file path, default: './spdata.dat'\n")
        if not output:
            output = "./spdata.dat"
    else:
        output = paras.output
    figs = args.figs if args.figs != "none" else FigType.none
    print("Reading gfile...")
    logging.info(f"Reading {filename}")
    data = read(filename)
    print("Begin mapping...")
    map_data, check_data = mapping_core(data, lsp, lst, psimax, figs)
    if figs in (FigType.save, FigType.show):
        check_plot(check_data, figs)
    print("Begin writing...")
    fortran_write(map_data, output)
    print(f"Finished writing to {output}")
