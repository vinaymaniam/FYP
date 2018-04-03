import numpy as np
import glob, os
import imageio
import cv2


def myglob(dir, pat):
    result = list()
    for file in os.listdir(dir):
        if file.endswith(pat):
            result.append(os.path.join(dir, file))
    return result


def load_images(paths):
    imgs = list()
    for i in range(0, len(paths), 1):
        x = imageio.imread(paths[i])
        if (len(x.shape) == 3):
            if (x.shape[2] == 3):
                x = rgb2ycbcr(x)
                x = x[:, :, 0]
        x = cv2.normalize(x.astype('float'), None, 0.0, 1.0, cv2.NORM_MINMAX)
        imgs.append(x)
    return imgs


def rgb2ycbcr(im_rgb):
    im_rgb = im_rgb.astype(np.float32)
    im_ycrcb = cv2.cvtColor(im_rgb, cv2.COLOR_RGB2YCR_CB)
    im_ycbcr = im_ycrcb[:, :, (0, 2, 1)].astype(np.float32)
    im_ycbcr[:, :, 0] = (im_ycbcr[:, :, 0] * (235 - 16) + 16) / 255.0  # to [16/255, 235/255]
    im_ycbcr[:, :, 1:] = (im_ycbcr[:, :, 1:] * (240 - 16) + 16) / 255.0  # to [16/255, 240/255]
    return im_ycbcr

# directory_x = 'Testing_Images/FRESH_upscaled/Set5'
# pattern = '.bmp'
# directory_y = 'Testing_Images/GT/Set5'
#
# XpathCell = myglob(directory_x, pattern )
# Xcell = load_images( XpathCell )
# YpathCell = myglob(directory_y, pattern )
# Ycell = load_images( YpathCell )
