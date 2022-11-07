import numpy as np

def coords_to_mask_2d (shape, coords, vals=True):
    """
    Create 2D mask from array of (x,y) coordinates
    :param shape:
    :param coords: Physical coordinates ([::-1] from heatmap coordinates)
    :param vals:
    """

    # Initialise and flatten mask
    mask_blank = np.zeros(shape)
    mask_line = mask_blank.reshape(-1)

    # Map coordinates to flattened heatmap
    c_line = coords[:, 0] * shape[1] + coords[:, 1]
    mask_line[c_line] = vals

    return mask_line.reshape(shape)
