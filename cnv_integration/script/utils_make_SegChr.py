#!/usr/bin/env python3

############################################################
# utils_make_SegChr.py
#
# Utility function to segment a chromosome into fixed windows.
############################################################

import math


def make_SegChr(CHR, start, end, window=1_000_000, PseudoChrNum=1_000_000_000):
    """
    Return:
      - Segmented_Chr: list[int]
      - Segmented_Chr_Len: list[float]
    """
    cont_start = (CHR * PseudoChrNum + start) / float(window)
    cont_end = (CHR * PseudoChrNum + end) / float(window)

    seg_start = math.floor(cont_start)
    seg_end = math.ceil(cont_end)

    # Fractional coverage of left and right boundary segments
    frac_left = math.ceil(cont_start) - cont_start
    frac_right = cont_end - math.floor(cont_end)

    segmented_chr = [i for i in range(seg_start, seg_end)]

    if int(cont_start) == int(cont_end):
        segmented_chr_len = [cont_end - cont_start]
    else:
        segmented_chr_len = (
            [frac_left]
            + [1 for _ in range(seg_start + 1, seg_end - 1)]
            + [frac_right]
        )

    return segmented_chr, segmented_chr_len

