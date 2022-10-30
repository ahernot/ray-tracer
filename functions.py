import numpy as np

def coords_to_mask_2d (shape, coords):
    mask_blank = np.zeros(shape, dtype=bool)
    mask_line = mask_blank.reshape(-1)

    c_line = coords[:, 0] * shape[0] + coords[:, 1]
    c_line_filtered = c_line[c_line < mask_line.shape[0]]

    mask_line[c_line_filtered] = True
    return mask_line.reshape(shape[::-1]).T
