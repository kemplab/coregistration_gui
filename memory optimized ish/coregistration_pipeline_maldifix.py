# import libraries

import numpy as np
import cv2


from scipy.optimize import dual_annealing
from scipy.spatial.distance import cdist


from skimage.filters import median, gaussian
from skimage.morphology import disk, label
from skimage.measure import regionprops_table
from skimage.feature import peak_local_max
from skimage.segmentation import watershed

from pyimzml.ImzMLParser import ImzMLParser

import pandas as pd
import os
import math
import time

import logging
import errno

from coregistration_utils import transform_image, mutual_information, crop, cleanup, relabel_random, \
    overlay_cells, visualize, get_connection_length, edge, ClusterImage, getmaxmin, record_reader, ASF, \
    getimagedimsfromcsv, readimagefromdf


"""
@author Tanya Roysam 
roysamt@gatech.edu
coregistration_pipeline.py
Coregistration pipeline for usage in GUI
"""



def preprocessing(maldi_path, confocal_path, pv, minsz = 2, maxsz = 5, bgremove = False):
    """
    Preprocessing steps of pipeline.
    Added morphological closing for usage with confocal image to improve homogeneity of confocal image

    :param maldi_path: (string) Path to MALDI MSI img
    :param confocal_path: (string) Path to confocal img
    :param pv: Shared-state variable ("progress var") for usage with progress bar.
    :param minsz: (int) min smallest dimension of cell in pixels
    :param maxsz: (int) max smallest dimension of cell in pixels
    :param bgremove: (bool) whether to remove background of confocal for alignment purposes only. Only recommended for
                    sparsely distributed cells, not clustered colonies.
    :return: (list) Preprocessed confocal, maldi images
    """
    f = 1

    # loading representative images
    logging.info('loaded imgs: '+ confocal_path + ' ' + maldi_path)
    logging.info('min cell size: '+str(minsz)+ 'px. max cell size: '+str(maxsz)+'px.')
    confocal = cv2.imread(confocal_path, -1)
    if len(confocal.shape) > 2:
        confocal = cv2.imread(confocal_path, 0)
        logging.warning("Multi-channel nuclear image was converted to grayscale. To preserve bit depth, use grayscale images only.")
    confocalbg = cv2.imread(confocal_path, 0)
    pv.set(pv.get()+f)
    maldi = cv2.imread(maldi_path, 0)
    pv.set(pv.get()+f)

    if maxsz < minsz:
        ms = maxsz
        maxsz = minsz
        minsz = ms
        logging.warning('Size upper limit was smaller than lower limit, limits were swapped')
    try:
        if confocal==None or confocalbg==None:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), confocal_path)
    except:
        pass
    try:
        if maldi==None:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), maldi_path)
    except:
        pass

    # normalizing image intensities
    confocal2 = cv2.normalize(confocal, None, alpha=0, beta=1, norm_type=cv2.NORM_MINMAX, dtype=cv2.CV_32F)
    pv.set(pv.get() + f)
    maldi2 = cv2.normalize(maldi, None, alpha=0, beta=1, norm_type=cv2.NORM_MINMAX, dtype=cv2.CV_32F)
    pv.set(pv.get() + f)

    # applying blurring filters to remove high-frequency noise and bring out the colony outlines
    # confocal2 = cv2.GaussianBlur(confocal2,(5,5),1)

    # Morphological Closing or Background Removal
    halfcs = maxsz #int(maxsz/2)
    if halfcs % 2 == 0: halfcs += 1

    if not bgremove:
        logging.info('closing with kernel size ('+str(halfcs)+','+str(halfcs)+') px, blurring')
        confocal2_closed = cv2.morphologyEx(confocal2, cv2.MORPH_CLOSE,
                                        cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (halfcs, halfcs)))
        pv.set(pv.get() + f)

        cinter = median(confocal2_closed, disk(5))
        confocal2 = gaussian(cinter, 1)
        pv.set(pv.get()+f)
    else:
        logging.info('removing background via normalizing, ASF from '+str(minsz)+' to '+str(maxsz)+' px, blur and thresh')
        norm_for_otsu = cv2.normalize(confocalbg, None, alpha=0, beta=255, norm_type=cv2.NORM_MINMAX)
        confocal2_asf = ASF(norm_for_otsu, [minsz, maxsz])

        confocal2_asf = cv2.GaussianBlur(confocal2_asf, (3, 3), 0)
        # Using Otsu's method to filter out the background
        mask_otsu = (1 - cv2.threshold(confocal2_asf, 0, 1, cv2.THRESH_OTSU)[1]).astype(bool)

        # Threshold gives binary image so using it as a mask over original image
        confocal2_thresh = np.copy(confocal2)
        confocal2_thresh[mask_otsu] = 0
        pv.set(pv.get() + f)

        confocal2 = confocal2_thresh
        pv.set(pv.get() + f)


    pv.set(pv.get()+f)
    # maldi2 = cv2.GaussianBlur(maldi, (5, 5), 1)
    maldi2 = gaussian(maldi2, 1)
    pv.set(pv.get()+f)


    # Output for GUI, Pipeline
    return [confocal2, maldi2]


def initial_alignment(confocal2, maldi2, xopt, pv, maldi_pixel_size_um, confocal_pixel_size_um):

    """
    Initial alignment of images

    :param confocal2: (ndarray, float) Preprocessed confocal image
    :param maldi2: (ndarray, float) Preprocessed MALDI MSI image
    :param xopt: (list, float) Initial aligment parameters
    :param pv: Shared-state variable ("progress var") for usage with progress bar.
    :param maldi_pixel_size_um: (float) Scale of MALDI MSI image in um/px
    :param confocal_pixel_size_um: (float) Scale of confocal image in um/px
    :return: (list, obj) initial alignment image, initial alignment score, transformed confocal & maldi, height, width,
            alignment parameters
    """
    # ORIGINAL

    # initial image dimensions
    h_conf, w_conf = confocal2.shape
    h_maldi, w_maldi = maldi2.shape
    logging.info('maldi scale: '+str(maldi_pixel_size_um)+' um/px. confocal scale: '+str(confocal_pixel_size_um)+' um/px.')

    if maldi_pixel_size_um < confocal_pixel_size_um:
        logging.warning('MALDI resolution higher than confocal. Check inputs if this is incorrect')



    m = h_conf / h_maldi
    # resizing confocal to match MALDI height
    confocal2 = cv2.resize(confocal2, dsize=(int(w_conf / m), h_maldi))
    h_conf_res, w_conf_res = confocal2.shape
    # matching the maldi and confocal width by leaving a stripe of zeros on the right for the narrower image
    h, w = max(h_conf_res, h_maldi), max(w_conf_res, w_maldi)  # common dimentions after shape matching
    logging.info('transforming maldi img')
    if w_conf_res > w_maldi:
        maldi3 = np.zeros((h, w))
        maldi3[:, :w_maldi] = maldi2
        confocal3 = confocal2
    else:
        confocal3 = np.zeros((h, w))
        confocal3[:, :w_conf_res] = confocal2
        maldi3 = maldi2
    if xopt[3] <= 0 or xopt[3]==1:
        xopt[3] = (maldi_pixel_size_um * h_maldi) / (confocal_pixel_size_um * h_conf)
    logging.info('initial alignment params, scaled for res: '+str(xopt))


    # setting boundaries for alignment parameters
    # the algorithm will look for a solution in parameter space bound by parameter+- parameter range
    # take a look at initial alignment, you may adjust initial alignment parameters in the box above
    im1, im2 = np.zeros((h, w, 3)), np.zeros((h, w, 3))
    im1[:, :, 0] = transform_image(maldi3, xopt[0], xopt[1], xopt[2], xopt[3])
    im1[:, :, 1] = transform_image(maldi3, xopt[0], xopt[1], xopt[2], xopt[3])
    im2[:, :, 2] = confocal3
    added_image = cv2.addWeighted(im1, 1, 2 * im2, 1, 0)
    sc = -round(mutual_information(xopt, confocal3, maldi3, pv=pv), 5)



    return [added_image, sc, confocal3, maldi3, h, w, xopt]

def final_alignment(confocal3, maldi3, h, w,
          project_dir, maldi_path, confocal_path,
          xopt = [0,0,0,1], maxiter = 1000, angle_range = 0.01,
          x_range = 100, y_range = 100,scale_range = 0.01, pv=None, csv=None, ions=[]):
    """
    Completes alignment and outputs aligned images in files specified.

    :param confocal3: (ndarray, float) Transformed confocal image
    :param maldi3: (ndarray, float) Transformed MALDI MSI image
    :param h: (int) height of both images
    :param w: (int) width of both images
    :param project_dir: (string) path to project directory
    :param maldi_path: (string) Path to original MALDI MSI image
    :param confocal_path: (string) Path to original confocal image
    :param xopt: (list, float) Alignment parameters
    :param maxiter: (int) Maximum iterations of dual annealing
    :param angle_range: (float) Angle range in radians within which to search for alignment
    :param x_range: (int) X range in px within which to search for alignment
    :param y_range: (int) Y range in px within which to search for alignment
    :param scale_range: (float) scale range within which to search for alignment
    :param pv: Shared-state variable ("progress var") for usage with progress bar.
    :param csv: (string) path to csv of raw MALDI data
    :param ions: (list, Float) list of ions
    :return: (list, obj) final aligned image, final alignment score, string describing results, (list, string) MALDI,
            confocal output paths, cropped confocal image, (list, ndarray) resized and cropped maldi images, (list,
            float) crop parameters for confocal images
    """

    # checking for problematic ranges & correcting
    if angle_range <= 0:
        logging.warning('Invalid angle range ' + str(angle_range) + ', angle range set to .01 rad')
        angle_range = .01
    if x_range <= 0:
        logging.warning('Invalid x range ' + str(x_range) + ', x range set to 1 px')
        x_range = 1
    if y_range <= 0:
        logging.warning('Invalid y range '+ str(y_range) + ', y range set to 1 px')
        y_range = 1
    if scale_range <= 0:
        logging.warning('Invalid scale range ' + str(scale_range) + ', scale range set to .01')



    logging.info('beginning dual annealment')
    p_opt = dual_annealing(mutual_information,
                           [(xopt[0] - angle_range, xopt[0] + angle_range), (xopt[1] - x_range, xopt[1] + x_range),
                            (xopt[2] - y_range, xopt[2] + y_range), (xopt[3] - scale_range, xopt[3] + scale_range)],
                           args=(confocal3, maldi3,pv), maxiter=maxiter, initial_temp=5e4, x0=xopt)

    finalsc = -round(p_opt.fun, 5)
    outstr = "Rotation (rad): " + str(round(p_opt.x[0], 5)) + " X shift: " + str(round(p_opt.x[1], 5)) + " Y shift: " + \
             str(round(p_opt.x[2], 5)) + " Scale: " + str(round(p_opt.x[3], 5))
    logging.info(outstr)
    # take a look at the alignment result
    maldi_transf = transform_image(maldi3, p_opt.x[0], p_opt.x[1], p_opt.x[2], p_opt.x[3])

    # cropping images for output
    [top_maldi, bottom_maldi, left_maldi, right_maldi] = crop(maldi_transf)
    [top_confocal, bottom_confocal, left_confocal, right_confocal] = crop(confocal3)
    top, bottom = max(top_maldi, top_confocal), min(bottom_maldi, bottom_confocal)
    left, right = max(left_maldi, left_confocal), min(right_maldi, right_confocal)
    im1, im2 = np.zeros((h, w, 3)), np.zeros((h, w, 3))
    im1[:, :, 0] = maldi_transf
    im1[:, :, 1] = maldi_transf
    im2[:, :, 2] = confocal3
    added_image = cv2.addWeighted(im1, 0.75, im2, 2, 0)
    finimg = added_image[top:bottom, left:right]

    confocal = cv2.imread(confocal_path, -1)
    if len(confocal.shape) > 2:
        confocal = cv2.imread(confocal_path, 0)
        logging.warning("Multi-channel nuclear image was converted to grayscale. To preserve bit depth, use grayscale images only.")

    im = cv2.imread(maldi_path, cv2.IMREAD_GRAYSCALE)

    h_conf, w_conf = confocal.shape
    h_maldi, w_maldi = im.shape
    m = h_conf / h_maldi

    h_conf_res, w_conf_res = confocal3.shape


    cropped_confocal = confocal[int(m * top):int(m * bottom), int(m * left):int(m * right)]

    confocal_crop_params = [m, top, bottom, left, right]

    # output aligned MALDI scaled up to original confocal resolution
    outdirs, resized_maldis = output_resize(project_dir, csv, ions, p_opt, cropped_confocal, confocal_crop_params,
                                            confocal_path, h_conf_res, h_maldi, w_conf_res, w_maldi)
    ydim, xdim = getimagedimsfromcsv(csv)
    resized_maldis = (csv, ydim, xdim, w_conf_res, w_maldi, h, w, xopt)
    return [finimg, finalsc, outstr, outdirs, cropped_confocal, resized_maldis, confocal_crop_params]


def output_resize(project_dir, csv, ions, p_opt, cropped_confocal, confocal_crop_params, confocal_path, h_conf_res,
                   h_maldi, w_conf_res, w_maldi):
    """
    Output post-alignment overlays for other confocal channels and MALDI ions

    :param project_dir: (string) directory of project
    :param csv: (string) path to csv of raw MALDI data
    :param ions: (list, float) list of ions
    :param p_opt: (OptimizeResult) object containing dual-annealed parameters
    :param cropped_confocal: (ndarray, float) cropped confocal image
    :param confocal_crop_params: (list, float) scale factor, top, bottom, left, right
    :param confocal_path: (string) path to DAPI image. other confocal channels should be found in the same directory
    :param h_conf_res: (int) height of resized confocal image
    :param h_maldi: (int) height of MALDI image
    :param w_conf_res: (int) width of resized confocal image
    :param w_maldi: (int) width of MALDI image
    :return: list: (list, string) output paths for GUI display, (list, ndarray) resized MALDI images
    """
    confocal = cv2.imread(confocal_path, -1)
    if len(confocal.shape) > 2:
        confocal = cv2.imread(confocal_path, 0)
        logging.warning("Multi-channel nuclear image was converted to grayscale. To preserve bit depth, use grayscale images only.")

    m, top, bottom, left, right = confocal_crop_params
    # output aligned MALDI scaled up to original confocal resolution

    ydim, xdim = getimagedimsfromcsv(csv)

    h, w = max(h_conf_res, h_maldi), max(w_conf_res, w_maldi)

    resized_maldis = []

    # create folders for output so that prior output is not overwritten

    maldi_output_path = project_dir + '/resized MALDI/'

    cn = 0
    while os.path.exists(maldi_output_path):
        cn += 1
        maldi_output_path = project_dir + "/resized MALDI ({})/".format(cn)
    os.mkdir(maldi_output_path)

    confocal_output = project_dir + '/resized confocal/'
    cn = 0
    while os.path.exists(confocal_output):
        cn += 1
        confocal_output = project_dir + "/resized confocal ({})/".format(cn)
    os.mkdir(confocal_output)

    # resize MALDI data
    maldi_df = pd.read_csv(csv)

    for ion in ions:
        name = "mz " + str(ion)

        try:
            im = readimagefromdf(ion, maldi_df, ydim, xdim)
        except KeyError:
            logging.warning('Ion {} not found in raw MALDI csv.'.format(ion))
            continue
        if w_conf_res > w_maldi:
            im2 = np.zeros((h, w))
            im2[:, :w_maldi] = im
        else:
            im2 = im
        im2 = transform_image(im2, p_opt.x[0], p_opt.x[1], p_opt.x[2], p_opt.x[3])
        im2 = im2[top:bottom, left:right]
        im2 = cv2.resize(im2, dsize=(cropped_confocal.shape[1], cropped_confocal.shape[0]))
        #resized_maldis.append(im2)
        cv2.imwrite(maldi_output_path + name + ".png", im2)

    del maldi_df # lorge

    conf_files = os.listdir(os.path.dirname(confocal_path))
    allowed_formats = [".png", ".tif", ".jpg", ".bmp"]
    m, top, bottom, left, right = confocal_crop_params

    # crop confocal & output

    for fname in conf_files:
        if fname[-4:] in allowed_formats:
            c = cv2.imread(os.path.dirname(confocal_path) + '\\' + fname, -1)
            if len(c.shape) > 2:
                c = cv2.imread(os.path.dirname(confocal_path) + '\\' + fname, 0)
                logging.warning(
                    "Multi-channel nuclear image was converted to grayscale. To preserve bit depth, use single-channel grayscale images only.")

            try:
                if (c.shape[0] >= confocal.shape[0]) & (c.shape[1] >= confocal.shape[1]):
                    cropped_c = c[int(m * top):int(bottom * m), int(m * left):int(right * m)]
                    cv2.imwrite(confocal_output + fname[:-4] + '_cropped.tif', cropped_c)
            except AttributeError:
                logging.warning('File ' + fname + ' was not readable.')

    logging.info('maldi output to ' + maldi_output_path)
    logging.info('confocal output to ' + confocal_output)

    #return [maldi_output_path, confocal_output], resized_maldis
    return [maldi_output_path, confocal_output], None


def pre_segmentation(maldi_path, conf_path, xopt, csv, low, high, project_dir, scattered_cells, remove_seeds,
            num_thresh, stds, ions, confocal_path, pv):
    """
    Preparing for segmentation in runs when existing alignment parameters are present

    :param maldi_path: (string) path to MALDI image
    :param conf_path: (string) path to confocal image
    :param xopt: (list, float) alignment parameters
    :param csv: (string) path to csv of MALDI values
    :param low: (int) Low range of radius in px
    :param high: (int) High range of radius in px
    :param project_dir: (string) path to project output
    :param scattered_cells: (bool) colony type
    :param remove_seeds: (bool) remove seeds that don't fall in colony mask
    :param num_thresh: (int) number of adaptive threshold window sizes
    :param stds: (float) number of standard deviations above which high intensity noise should be thresholded
    :param ions: (list, float) MALDI m/z values to overlay
    :param confocal_path: (string) path to confocal image
    :param pv: shared-state progress bar variable
    :return: list: (ndarray) labeled nuclei, (DataFrame) cell labels/metrics, cropped confocal image
    """

    confocal = cv2.imread(conf_path, -1)
    if len(confocal.shape) > 2:
        confocal = cv2.imread(confocal_path, 0)
        logging.warning("Multi-channel nuclear image was converted to grayscale. To preserve bit depth, use single-channel grayscale images only.")
    im = cv2.imread(maldi_path, cv2.IMREAD_GRAYSCALE)

    ydim, xdim = getimagedimsfromcsv(csv)


    h_conf, w_conf = confocal.shape
    h_maldi, w_maldi = im.shape
    m = h_conf / h_maldi


    confocal2 = cv2.resize(confocal, dsize=(int(w_conf / m), h_maldi))
    h_conf_res, w_conf_res = confocal2.shape
    # matching the maldi and confocal width by leaving a stripe of zeros on the right for the narrower image
    h, w = max(h_conf_res, h_maldi), max(w_conf_res, w_maldi)  # common dimentions after shape matching
    logging.info('transforming maldi img')
    if w_conf_res > w_maldi:
        maldi3 = np.zeros((h, w))
        maldi3[:, :w_maldi] = im
        confocal3 = confocal2
    else:
        confocal3 = np.zeros((h, w))
        confocal3[:, :w_conf_res] = confocal2
        maldi3 = im


    maldi_transf = transform_image(maldi3, xopt[0], xopt[1], xopt[2], xopt[3])
    # cropping images for output
    [top_maldi, bottom_maldi, left_maldi, right_maldi] = crop(maldi_transf)
    [top_confocal, bottom_confocal, left_confocal, right_confocal] = crop(confocal3)

    top, bottom = max(top_maldi, top_confocal), min(bottom_maldi, bottom_confocal)
    left, right = max(left_maldi, left_confocal), min(right_maldi, right_confocal)



    cropped_confocal = confocal[int(m * top):int(m * bottom), int(m * left):int(m * right)]

    confocal_crop_params = [m, top, bottom, left, right]
    resized_maldis = []
    # output aligned MALDI scaled up to original confocal resolution

    resized_maldis = (csv, ydim, xdim, w_conf_res, w_maldi, h, w, xopt)
    """for ion in ions:
        try:
            im = readimagefromdf(ion, maldi_df, ydim, xdim)
        except KeyError:
            logging.warning('Ion {} not found in raw MALDI csv.'.format(ion))
            continue
        if w_conf_res > w_maldi:
            im2 = np.zeros((h, w))
            im2[:, :w_maldi] = im
        else:
            im2 = im
        im2 = transform_image(im2, xopt[0], xopt[1], xopt[2], xopt[3])
        im2 = im2[top:bottom, left:right]
        im2 = cv2.resize(im2, dsize=(cropped_confocal.shape[1], cropped_confocal.shape[0]))
        resized_maldis.append(im2)"""


    return segment(cropped_confocal, low, high, project_dir, scattered_cells, remove_seeds, num_thresh, stds,
                   resized_maldis, ions, confocal_crop_params, confocal_path, pv)


def output_segmentation(cropped_confocal, resized_maldis, cells, metrics_output_folder, ions, csv_output, confocal_crop_params, confocal_path, pv):
    """
    Output/overlaying segmented data
    :param cropped_confocal: (ndarray) cropped confocal image
    :param resized_maldis: (list, ndarray) resized maldi images
    :param cells: (DataFrame) connected components
    :param metrics_output_folder: (string) path to output metrics in
    :param confocal_output: (string) path to output cropped confocal images
    :param ions: (list, float) ions
    :param csv_output: (string) path to output csv of prelim metrics
    :param confocal_crop_params: (list, float) crop parameters for confocal images
    :param confocal_path: (string) path to original confocal image
    :param pv: Shared-state variable ("progress var") for usage with progress bar.
    :return: (list, Figure) list of confocal overlay figures, (list, Figure) list of maldi overlay figures
    """
    # PV is = 2 + len(ions)
    f = 1
    Y = np.array(cells['centroid-1'])
    X = np.array(cells['centroid-0'])
    area = np.array(cells['area'])
    perimeter = np.array(cells['perimeter'])
    eccentricity = np.array(cells['eccentricity'])

    # overlay the dapi channel
    dapi_overlaid = overlay_cells(cropped_confocal, cells,
                                    'ASF.png')[0]

    excel_header = ["X", "Y", "Area", "Perimeter", "Eccentricity", "dapi"]
    excel_data = [Y, np.max(X) - X, area, perimeter, eccentricity, dapi_overlaid]
    conf_figs = []
    maldi_figs = []

    cofn = metrics_output_folder+"/confocal overlays/"
    cn = 0
    while os.path.exists(cofn):
        cn += 1
        cofn = metrics_output_folder+"/confocal overlays ({})/".format(cn)
    os.mkdir(cofn)

    original_confocal = cv2.imread(confocal_path, -1)
    if len(original_confocal.shape) > 2:
        original_confocal = cv2.imread(confocal_path, 0)
        logging.warning("Multi-channel nuclear image was converted to grayscale. To preserve bit depth, use grayscale images only.")

    conf_files = os.listdir(os.path.dirname(confocal_path))
    allowed_formats = [".png",".tif",".jpg",".bmp"]
    m, top, bottom, left, right = confocal_crop_params

    for fname in conf_files:
        if fname[-4:] in allowed_formats:
            confocal = cv2.imread(os.path.dirname(confocal_path)+'/'+fname,-1)
            if len(confocal.shape) > 2:
                confocal = cv2.imread(os.path.dirname(confocal_path)+'/'+fname, 0)
                logging.warning(
                    "Multi-channel nuclear image was converted to grayscale. To preserve bit depth, use grayscale images only.")
            if (confocal.shape[0] >= original_confocal.shape[0]) & (confocal.shape[1] >= original_confocal.shape[1]):
                cropped_confocal = confocal[int(m * top):int(bottom * m), int(m * left):int(right * m)]
                try:
                    overlay_conf_out = overlay_cells(cropped_confocal, cells,
                                                 cofn + fname[:-4] + ".png")
                    excel_data.append(overlay_conf_out[0])
                    excel_header.append(fname[:-4])
                    #conf_figs.append(overlay_conf_out[2])
                    logging.info(fname + ' cropped and overlaid')
                except ValueError:
                    logging.warning('non-grayscale image in confocal folder (cannot overlay): ' + fname)

            else:
                logging.warning('File '+ fname+ ' found in confocal input directory that was not compatible')

    pv.set(pv.get() + f)

    mn = 0
    moutpath = metrics_output_folder+"/MALDI overlays/"
    while os.path.exists(moutpath):
        mn += 1
        moutpath =  metrics_output_folder+"/MALDI overlays ({})/".format(mn)
    os.mkdir(moutpath)

    #resized_maldis = []
    # output aligned MALDI scaled up to original confocal resolution


    # output aligned MALDI scaled up to original confocal resolution
    csv, ydim, xdim, w_conf_res, w_maldi, h, w, xopt = resized_maldis
    maldi_df = pd.read_csv(csv)

    for ind, ion in enumerate(ions):
        try:
            im = readimagefromdf(ion, maldi_df, ydim, xdim)
        except KeyError:
            logging.warning('Ion {} not found in raw MALDI csv.'.format(ion))
            continue
        if w_conf_res > w_maldi:
            im2 = np.zeros((h, w))
            im2[:, :w_maldi] = im
        else:
            im2 = im
        im2 = transform_image(im2, xopt[0], xopt[1], xopt[2], xopt[3])
        im2 = im2[top:bottom, left:right]
        im2 = cv2.resize(im2, dsize=(cropped_confocal.shape[1], cropped_confocal.shape[0]))

        name = "mz " + str(ions[ind])
        overlay_out = overlay_cells(im2, cells, moutpath + name + ".png")
        pv.set(pv.get() + f)
        excel_data.append(overlay_out[0])
        excel_header.append(name)
        # maldi_figs.append(overlay_out[2])
        logging.info(str(ions[ind]) + ' cropped and overlaid')

        #resized_maldis.append(im2)

    del maldi_df

    """for ind, im in enumerate(resized_maldis):
        # same order as ions.
        name = "mz " + str(ions[ind])
        overlay_out = overlay_cells(im, cells, moutpath+name+".png")
        pv.set(pv.get() + f)
        excel_data.append(overlay_out[0])
        excel_header.append(name)
        #maldi_figs.append(overlay_out[2])
        logging.info(str(ions[ind]) + ' cropped and overlaid')
    """


    out_dict = {}

    for n, header in enumerate(excel_header):
        out_dict.update({header:excel_data[n]})

    out = pd.DataFrame(out_dict)
    out.to_csv(csv_output)
    logging.info('labels output to '+csv_output)
    pv.set(pv.get() + f)


    #return conf_figs, maldi_figs

def LoG(low, high, step, confocal_2, pv):
    """
    Laplacian of Gaussian implementation.
    'Improved Automatic Detection and Segmentation of Cell Nuclei in Histopathology Images', Al - Kofahi
    et al(2010) https://ieeexplore.ieee.org/abstract/document/5306149

    :param low: (int) low end of radius range
    :param high: (int) high end of radius range
    :param step: (int) determines how many radii within range are tested
    :param confocal_2: (ndarray) confocal image on which LoG is computed
    :param pv: shared-state variable for progress bar
    :return: max_LoG: (ndarray) same size as confocal image, contains image used for seeding
    """
    # calculate range of sigmas based on given radii
    # sigma = radius / sqrt(2)
    sigma = np.arange(low, high + step, step) / math.sqrt(2)

    logging.info('beginning laplacian of gaussian, sigmas ' + str(sigma))

    #all_LoG = np.zeros((confocal_2.shape[0], confocal_2.shape[1], len(sigma)))
    max_LoG = np.zeros((confocal_2.shape[0], confocal_2.shape[1]))

    # laplacian of gaussian
    for i, s in enumerate(sigma):

        im = gaussian(confocal_2, s)
        im = cv2.Laplacian(im, cv2.CV_64F, ksize=1)

        # Normalizing by multiplying by squared sigma
        n_im_norm = im * (s ** 2)
        if i == 0:
            max_LoG = n_im_norm

        # Adding to array storing all LoG results
        #all_LoG[:, :, i] = n_im_norm
        max_LoG = np.amin([n_im_norm, max_LoG], axis=0)

        pv.set(pv.get() + 1)

    # Finding minimum value of each pixel across all LoG results at different sigma values
    #max_LoG = np.amin(all_LoG, axis=2)
    # Inverting so blob centers are local maxima rather than minima
    max_LoG = np.max(max_LoG) - max_LoG

    return max_LoG

def cell_mask(confocal_thresh, confocal_original, high, num_thresh, pv):
    """
    Generate mask of all cells using multiscale adaptive threshold

    :param confocal_thresh: (ndarray) confocal image to threshold
    :param confocal_original: (ndarray) original confocal image
    :param high: (int) in px, highest possible size of nucleus
    :param num_thresh: (int) number of thresholds. higher number -> more foreground
    :param pv: shared-state variable for progress bar
    :return: total_mask: (ndarray) mask of all the nuclei
    """

    # Gaussian blurring in preparation for adaptive thresh
    confocal_blur = cv2.normalize(cv2.GaussianBlur(confocal_thresh, (3, 3), 0), None, alpha=0, beta=255,
                                  norm_type=cv2.NORM_MINMAX, dtype=cv2.CV_8UC1)
    pv.set(pv.get() + 1)

    total_mask = np.zeros(confocal_original.shape)

    # Adaptive threshold at multiple window sizes
    # Minimum window size is diameter of largest cell
    # And maximum window size is the minimum shape of the confocal image
    # Choosing 10 equally spaced kernel sizes between these two
    # Increase this number if too few cells are registered
    logging.info('beginning thresholding')
    for sz in np.linspace(high * 2, np.min(confocal_original.shape), num_thresh):
        # Rounding floats and making sure they are odd
        sz = int(round(sz))
        if sz % 2 == 0: sz += 1
        confocal_mask = cv2.adaptiveThreshold(confocal_blur, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, sz,
                                              0)
        pv.set(pv.get() + 1)
        # Bitwise or to combine this mask with the total mask
        total_mask = np.bitwise_or(total_mask.astype(np.uint8), confocal_mask.astype(np.uint8))
        pv.set(pv.get() + 1)
    #cv2.imwrite("C:\\Users\\tanya\\PycharmProjects\\gui\\coregistration-master\\for github\\mask.tif", total_mask)
    return total_mask

def seed(max_LoG, min_dist, watershed_img, total_mask, remove_seeds=False, m=None, pv=None):
    """
    seed and segment image according to LoG results

    :param max_LoG: (ndarray) output of LoG, maximum value of LoG across all scales at each pixel
    :param min_dist: (int) minimum distance between two seeds
    :param watershed_img: (ndarray) image to watershed
    :param total_mask: (ndarray) mask image to watershed
    :param remove_seeds: (bool) remove seeds outside colony
    :param m: (ndarray) mask of colony
    :param pv: shared-state progress variable
    :return: nuclei_labels: (ndarray) labeled image of nuclei
    """
    # Finding local maxima. Minimum distance between maxima is minimum radius of nucleus
    logging.info('finding seeds')
    loc_inds = peak_local_max(max_LoG, min_distance=min_dist)
    pv.set(pv.get() + 1)

    locs = np.zeros(max_LoG.shape)
    locs[loc_inds[:, 0], loc_inds[:, 1]] = 1
    # From before: m is an approximate mask of colony. This removes false seeds generated by noise
    if remove_seeds:
        locs[m == 0] = 0

    # Label and randomize label of seeds
    labels = label(locs)
    pv.set(pv.get() + 1)
    labels_rand = relabel_random(labels)
    pv.set(pv.get() + 1)

    # Watershed
    logging.info('watershedding')
    nuclei_labels = watershed(watershed_img, labels_rand, mask=total_mask)
    pv.set(pv.get() + 1)

    return nuclei_labels

def segment_with_LoG(cropped_confocal, high, low, stds, num_thresh, remove_seeds, pv):
    """
    Segment image using LoG. for clustered colonies

    :param cropped_confocal: (ndarray) cropped DAPI image to segment
    :param high: (int) high range of radii
    :param low: (int) low range of radii
    :param stds: (float) standard deviations above which to threshold high-intensity noise
    :param num_thresh: (int) number of thresholds for adaptive threshold. increase to increase fg
    :param remove_seeds: (bool) remove seeds outside colony. recommend True
    :param pv: shared-state progress bar variable
    :return: ASF_conf: (ndarray) post-ASF confocal, nuclei_labels: (ndarray) labeled image of nuclei
    """
    logging.info('segmenting as clustered colony')
    # FIXME change this for memory/speed? default is 1, increase to reduce alloc
    step = 2

    conf = cropped_confocal

    gbdim = high * 6 + 1

    # For dividing
    clow = int(low / 2)
    chigh = int(low)

    # Finding the median and standard deviation of intensities of the confocal image
    med = np.median(conf)
    std = np.std(conf)

    # Thresholding outliers out
    # High-intensity outliers distort results so I truncate them
    thr = int(round(med + stds * std))
    lowthr = int(round(med - stds * std))
    if lowthr < 0: lowthr = 0

    confocal_2 = cv2.threshold(conf, thr, 255, cv2.THRESH_TRUNC)[1]

    confocal_3 = cv2.GaussianBlur(conf, (gbdim, gbdim), 0)
    confocal_3 = cv2.subtract(conf, confocal_3)

    ASF_conf = ASF(confocal_3, range(clow, chigh, 1))

    pv.set(pv.get() + 1)

    # Generating approximate mask of entire colony
    # Useful for later, deciding which seeds are cells and which are noise
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (high * 4 + 1, high * 2 + 1))
    closed = cv2.morphologyEx(confocal_2, cv2.MORPH_CLOSE, kernel, iterations=1)
    pv.set(pv.get() + 1)
    m = cv2.threshold(closed, 0, 255, cv2.THRESH_OTSU + cv2.THRESH_BINARY)[1]
    m = cv2.morphologyEx(m, cv2.MORPH_OPEN, kernel, iterations=1)

    pv.set(pv.get() + 1)

    # find mask of colony
    total_mask = cell_mask(confocal_3, conf, high, num_thresh, pv)

    # As in paper- standard deviation sigma of gaussian = radius/sqrt(2)
    # Finding range of sigmas using range of radii
    max_LoG = LoG(low, high, step, confocal_2, pv)

    nuclei_labels = seed(max_LoG, low, total_mask, total_mask, remove_seeds, m, pv)

    return ASF_conf, nuclei_labels

def segment_with_ASF_LoG(cropped_confocal, low, high, pv):
    """
    Segment image using Alternating Sequential Filter, LoG. for scattered colonies

    :param cropped_confocal: (ndarray) cropped confocal image
    :param low: (int) low range of radii
    :param high: (int) high range of radii
    :param pv: shared-state progress bar variable
    :return: ASF_conf: (ndarray) post-ASF confocal, nuclei_labels: (ndarray) labeled image of nuclei
    """
    conf = cropped_confocal
    step = 1
    norm_for_otsu = cv2.normalize(conf, None, alpha=0, beta=255, norm_type=cv2.NORM_MINMAX)
    pv.set(pv.get() + 1)
    confocal2_asf = ASF(norm_for_otsu, [low * 2, high * 2])
    pv.set(pv.get() + 1)

    confocal2_asf = cv2.GaussianBlur(confocal2_asf, (3, 3), 0)
    ASF_conf = confocal2_asf
    pv.set(pv.get() + 1)
    # Using Otsu's method to filter out the background
    retval, otsu = cv2.threshold(confocal2_asf, 0, 1, cv2.THRESH_OTSU)
    mask_otsu = (1 - otsu).astype(bool)
    pv.set(pv.get() + 1)

    logging.info('ASF, Otsu\'s thresh at retval ' + str(retval))

    # Threshold gives binary image so using it as a mask over original image
    confocal2_thresh = np.copy(conf)
    confocal2_thresh[mask_otsu] = 0

    # As in paper- standard deviation sigma of gaussian = radius/sqrt(2)
    # Finding range of sigmas using range of radii
    max_LoG = LoG(low, high, step, confocal2_thresh, pv)
    pv.set(pv.get() + 1)

    nuclei_labels = seed(max_LoG, high, confocal2_thresh, otsu, pv=pv)

    return ASF_conf, nuclei_labels


def segment(cropped_confocal, low, high, project_dir, scattered_cells, remove_seeds, num_thresh, stds, resized_maldis,
          ions, confocal_crop_params, confocal_path, pv):
    """
    Segment and overlay images
    :param cropped_confocal: (ndarray) cropped confocal image
    :param low: (int) Low range of radius in px
    :param high: (int) Upper range of radius in px
    :param project_dir: (string) path to project fold
    :param scattered_cells: (bool) if the cells are separated or clustered. determines preprocessing method
    :param remove_seeds: (bool) if the seeds outside the general colony shape should be removed, only matters if
    clustered cells
    :param num_thresh: (int) Number of adaptive thresholds. More thresholds -> more cell signal picked up
    :param stds: (float) Signal above specified standard deviations from mean are ceil
    :param resized_maldis: (list, ndarray) list of resized maldi images
    :param ions: (list, float) list of ions
    :param confocal_crop_params: (list, float) crop parameters for confocal images
    :param confocal_path: (string) path to original confocal images
    :param pv: Shared-state variable ("progress var") for usage with progress bar.
    :return: list: (ndarray) labeled nuclei, (DataFrame) cell labels/metrics, cropped confocal image
    """
    # If scattered_cells == False:
    # pv = 3 + 2*num_thresh + (high - low + 1) + 5 + 3 + 2 + len(ions)= 14 + 2*num_thresh + (high - low) + len(ions)
    # if true
    # pv = 3 + high - low + 1 + 5 + 3 + 2 + len(ions)= 14 + high - low + len(ions)
    f = 1
    # start cell segmentation



    logging.info('Beginning segmentation run.')
    logging.info('min radius: '+str(low)+'px. max radius: '+str(high)+'px.')
    logging.info('Num thresholds: '+str(num_thresh))
    logging.info('Upper threshold: '+str(stds)+ ' standard deviations above mean intensity')

    # correcting problematic inputs

    if low > high:
        lt = low
        low = high
        high = lt
        logging.warning('Lower radius limit was greater than upper radius limit, limits swapped')

    num_thresh = abs(num_thresh)
    if num_thresh <= 1:
        num_thresh = 2
        logging.warning('Number of thresholds <=1 was specified. Must have at least 2 thresholds.')

    if num_thresh%1 > 0:
        num_thresh = round(num_thresh)
        logging.warning('Number of thresholds must be positive integer >= 2.')


    """Laplacian of Gaussian method from:
    'Improved Automatic Detection and Segmentation of Cell Nuclei in Histopathology Images', Al - Kofahi
    et al(2010) https://ieeexplore.ieee.org/abstract/document/5306149"""
    if not scattered_cells:
        ASF_conf, nuclei_labels = segment_with_LoG(cropped_confocal, high, low, stds, num_thresh, remove_seeds, pv)
    else:
        ASF_conf, nuclei_labels = segment_with_ASF_LoG(cropped_confocal, low, high, pv)



    logging.info('cleaning labels by size, overlaying')
    cellprops = pd.DataFrame(regionprops_table(nuclei_labels,
                                               properties=('label','centroid','area','perimeter','eccentricity', 'coords')))
    pv.set(pv.get() + f)


    mn = 0
    seg_save = project_dir + '/segmentation out/'
    while os.path.exists(seg_save):
        mn += 1
        seg_save = project_dir + "/segmentation out ({})/".format(mn)
    os.mkdir(seg_save)

    lblout =  seg_save + 'labels.csv'


    logging.info('saving raw labels: '+lblout)
    cellprops.to_csv(lblout)
    pv.set(pv.get() + f)
    cells, X, Y, area, perimeter, eccentricity = cleanup(cellprops)
    pv.set(pv.get() + f)

    output_segmentation(ASF_conf, resized_maldis, cells, seg_save, ions, seg_save + 'metrics.csv', confocal_crop_params,
                        confocal_path, pv)

    points = np.concatenate(cells['coords'].tolist())
    nuclei_mask = np.zeros(nuclei_labels.shape)
    nuclei_mask[points[:,0], points[:,1]] = 1
    nuclei_labels[np.where(nuclei_mask==0)] = 0

    #return [nuclei_labels, cells, cropped_confocal, seg_save, confocal_crop_params, resized_maldis]
    return [nuclei_labels, cells, cropped_confocal, seg_save, confocal_crop_params, None]

def manual_km(coords, c1, c2):
    cd = cdist(coords, [c1, c2])
    return np.argmax(cd, axis=1)

def cell_network(cells_path, metrics_range, num_value_layers, min_clust_size, connection, project_dir, edge_margin,
          cropped_confocal, labels_df, vis, pv):
    """
    Calculate cell network metrics

    :param cells_path: (string) path to preliminary metrics
    :param metrics_range: (int) range of cells
    :param num_value_layers: (int) number of intensity classes for kmeans clustering
    :param min_clust_size: (int) minimum size of kmeans cluster
    :param connection: (float) length of connection between cells. if it is -1 it is calculated automatically
    :param project_dir: (string) path to project directory
    :param edge_margin: (int) px margin to exclude when checking for edge cells
    :param cropped_confocal: (ndarray) cropped confocal image
    :param labels_df: (DataFrame) df containing coordinate labels of each nucleus for visualizations
    :param vis: (bool) if True, will generate visualizations
    :param pv: Shared-state variable ("progress var") for usage with progress bar.
    :return: (list, string) Path of cell network metrics csv output
    """

    # loading in data
    prelim_metricspath = cells_path + '/' + 'metrics.csv'
    try:
        cells = pd.read_csv(prelim_metricspath)
    except FileNotFoundError:
        cnm_outpath = 'Preliminary metrics file generated in segmentation step was not found in directory.\n' \
                  'Please refrain from moving, modifying, or deleting \'metrics.csv\', found in the same\n' \
                  'directory as the label csv, until the application is closed.'

        raise FileNotFoundError(errno.ENOENT, cnm_outpath, prelim_metricspath)

    try:
        cells = cells.drop(columns='Unnamed: 0')
        variable_headers = cells.columns[5:]
    except KeyError:
        variable_headers = cells.columns[6:]

    num_cells = len(cells)
    num_variables = cells.shape[1] - 5

    # x and y coordinates of cells
    X = np.array(cells['X'])
    Y = np.array(cells['Y'])

    # connection < 0 indicates connection length should be calculated
    if connection <= 0:
        connection = get_connection_length(X, Y)  # standardize to 20 for testing
        logging.info('Generating cell-cell connection length: ' + str(connection))
    pv.set(pv.get()+1)
    region = connection * abs(int(round(metrics_range)))

    # Finding euclidean distance between each coordinate. If there are n coordinates, dists is nxn 2 dimensional array
    coords = np.array([X, Y]).T
    logging.info('finding distances between cells')
    dists = cdist(coords, coords, metric='euclidean')
    pv.set(pv.get() + 1)

    # Finding the global list of cell neighborhoods
    neighborhoods = dists < region
    neighborhoods[range(len(neighborhoods)), range(len(neighborhoods))] = 0
    connected = np.bitwise_and(neighborhoods, dists < connection)
    logging.info('finding neighborhood regions of cells')
    global_neigh_region = [np.where(x) for x in neighborhoods]
    pv.set(pv.get() + 1)

    # Finding the number of neighbors
    # find neighbors (within region) and connections (neighbors within connection distance)
    num_neighbors = np.sum(neighborhoods.astype(int), 0)
    num_connections = np.sum(connected.astype(int), 0)

    has_connections = num_connections > 0
    avgLength = np.nan_to_num(np.sum(connected * dists, 0) / num_connections)
    avgLength[has_connections != True] = connection
    density = np.nan_to_num(num_connections / avgLength)

    # distance from self to nearest neighbor
    nearest_neighbor_dist = np.ma.masked_equal(dists, 0.0, copy=False).min(axis=0)

    # global list of neighbors omitting those with no neighbors/no connections
    nz_mask = np.array(global_neigh_region, dtype=object)

    nonzero_connections = nz_mask[has_connections]
    nonzero_connections = nonzero_connections.tolist()

    cells['number of connections'] = num_connections
    cells['density'] = density
    cells['average length'] = avgLength
    cells['nearest_neighbor_dist'] = nearest_neighbor_dist

    averageNeigh = np.zeros((num_variables, num_cells))
    selfOverAvg = np.zeros((num_variables, num_cells))
    stdNeigh = np.zeros((num_variables, num_cells))
    rangeNeigh = np.zeros((num_variables, num_cells))
    avgCluster = np.zeros((num_variables, num_cells))
    clusterSize = np.zeros((num_variables, num_cells))
    selfOverAvgCluster = np.zeros((num_variables, num_cells))

    # 33 s for above, 2s/variable loop

    # for each variable, calculate absolute and neighbor-relative metrics
    for ind, var in enumerate(variable_headers):
        try:
            logging.info('calculating metrics for '+var)

            # Finding values
            values = np.array(cells[var])
            values_nz = values[has_connections]

            # values within neighborhood
            neigh_values = [values[tuple(locs)] for locs in nonzero_connections]

            # finding mean and std of neighborhoods of cells w/ connections
            mwi = np.array([np.mean(nval) for nval in neigh_values])
            stdwi = np.array([np.std(nval) for nval in neigh_values])

            # neighbor relative value of neighborhoods of cells w connections
            # "self over average"
            soa = values_nz / mwi
            mask1 = np.bitwise_and(values_nz == 0, mwi == 0)
            soa[mask1] = 1
            mask10 = np.bitwise_and(mwi == 0, values_nz != 0)
            soa[mask10] = 10

            # finding ranges of each neighborhood of cells w connections
            range_nz = np.array([np.ptp(values[tuple(locs)]) for locs in nonzero_connections])

            # mapping range, mean, std. and soa to arrays including cells with no connections
            ranges = np.zeros(num_cells)
            ranges[has_connections] = range_nz
            mean = np.array(values)
            mean[has_connections] = mwi
            std = np.zeros(num_cells)
            std[has_connections] = stdwi
            soaverage = np.ones(num_cells)
            soaverage[has_connections] = soa

            averageNeigh[ind] = mean
            selfOverAvg[ind] = soaverage
            stdNeigh[ind] = std
            rangeNeigh[ind] = ranges

            # storing in df
            cells[var + '_average'] = mean
            cells[var + '_std'] = std
            cells[var + '_self_value_over_average'] = soaverage
            cells[var + '_range'] = ranges
        except:
            logging.warning('Error encountered for {}-possible NaN/large value'.format(var))

    # around 50s for basic
    pv.set(pv.get() + 1)

    # calculate which cells are on the edge of a colony
    # an edge cell is defined by having at least one side free of neighbors
    # one "side" is defined as pi/2 rad
    logging.info('calculating which cells are edge cells')
    edge_vals = edge(X, Y, global_neigh_region, has_connections, edge_margin)
    cells['edge'] = edge_vals

    pv.set(pv.get() + 1)

    # Calculate whether cells are dividing or not
    # Three features: neighbor-relative DAPI intensity, area, and distance to nearest neighbor
    # All features are minmax normalized from 0 to 1
    dapisoa = np.array(cells['dapi_self_value_over_average'])
    dapisoa = (dapisoa - min(dapisoa)) / (max(dapisoa) - min(dapisoa))
    dapisoa *= 2 # weighting this the highest

    area = np.array(cells['Area'])
    area = (area - min(area)) / (max(area) - min(area))

    nearest_neighbor_dist = (nearest_neighbor_dist - min(nearest_neighbor_dist)) / \
                            (max(nearest_neighbor_dist) - min(nearest_neighbor_dist))
    # future: add inputs to change this
    # kmeans centers for dividing (dc) and not dividing (ndc) cells
    # order: dapi_SoA, area, nearest neighbor distance
    # these were determined using control images under no pi3k inhibition and validated using p-AKT stain
    dc = [0.92918492, 0.21520196, 0.18716013]
    ndc = [0.21033933, 0.24907662, 0.17722891]

    pts = np.array([dapisoa, area, nearest_neighbor_dist])

    # KMeans clustering to find which ones are dividing
    div_lbls = manual_km(pts.T, dc, ndc)
    cells['dividing'] = div_lbls


    # Clustering: 25 s, 2.78 s per variable
    for ind, var in enumerate(variable_headers):
        try:
            logging.info('clustering '+var)
            clusterSize[ind, :], avgCluster[ind, :] = ClusterImage(averageNeigh[ind, :], num_value_layers, X, Y, connection,
                                                                   min_clust_size)
            selfOverAvgCluster[ind, :] = averageNeigh[ind, :] / avgCluster[ind, :]

            cells[var + '_cluster_size'] = clusterSize[ind, :]
            cells[var + '_average_cluster'] = avgCluster[ind, :]
            cells[var + '_self_value_over_average_cluster'] = selfOverAvgCluster[ind, :]
        except ValueError:
            logging.warning("{} contains NaN, could not cluster".format(var))
    pv.set(pv.get() + 1)
    min_edge_dist = 100
    
    
    edgeDist, centerDist = np.zeros(num_cells), np.zeros(num_cells)

    edge_coords = np.array([X[edge_vals == 1], Y[edge_vals == 1]]).T

    logging.info('calculating edge dists')

    edgeDist[edge_vals == 1] = 0
    edge_dists = cdist(coords, edge_coords, metric='euclidean')  # dim 1 is coords, dim 2 is edge
    min_edge_dists = np.min(edge_dists, 1)
    min_edge_dists[min_edge_dists < min_edge_dist] = min_edge_dist
    edgeDist[edge_vals == 0] = min_edge_dists[edge_vals == 0]

    # calculate distance to the center
    max_center_dist = 1 / min_edge_dist ** 0.25

    centerDist[edge_vals == 1] = max_center_dist
    centerDist[edge_vals == 0] = 1 / edgeDist[edge_vals == 0] ** .25

    cells['edgeDist'] = edgeDist
    cells['centerDist'] = centerDist
    pv.set(pv.get() + 1)
    cnm_outpath = project_dir + '/cell_network_metrics.csv'

    cells.to_csv(cnm_outpath)

    # visualize dividing cells regardless of visualization preference
    med = np.median(cropped_confocal)
    std = np.std(cropped_confocal)

    thr = int(round(med + std))
    lowthr = int(round(med - std))
    if lowthr < 0: lowthr = 0

    dividing_mask = np.array(div_lbls) == 1

    div_lbls += max(div_lbls)

    visualize(div_lbls*3+1, labels_df, cropped_confocal.shape[0], cropped_confocal.shape[1], project_dir + '/dividing_labeled.png')

    # visualization takes a while so it is optional
    if vis:
        vispath = project_dir + '/visualizations/'
        vpn = 0
        while os.path.exists(vispath):
            vpn += 1
            vispath = project_dir + '/visualizations ({})/'.format(vpn)
        os.mkdir(vispath)
        logging.info('generating visualizations')
        t = time.time()
        cnm_out(cells, labels_df, cropped_confocal.shape[0], cropped_confocal.shape[1], vispath)
        elapsed = time.time()-t
        logging.info('visualizations took '+str(int(elapsed/60))+' m to generate')


    pv.set(pv.get() + 1)

    return [cnm_outpath]

def cnm_out(metrics_df, labels_df, h, w, vispath):
    """
    Calculate visualizations of cell network metrics
    :param metrics_df: (DataFrame) df with calculated metrics
    :param labels_df: (DataFrame) df with coordinate labels of nuclei
    :param h: (int) Height of img
    :param w: (int) Width of img
    :param vispath: (string) Path to store visualizations in (directory)
    :return: None
    """
    #7th col onwards is vis
    names = metrics_df.columns[6:]
    for n in names:
        vals = np.array(metrics_df[n])
        visualize(vals, labels_df, h, w, vispath+n+'.png')

def imzml_parse(MALDI_path, ions, project_dir, pv):
    """
    Parse ImzML file and generate ion images, raw ion data csv

    :param MALDI_path: (string) path to ImzML file
    :param ions: (list, float) List of ion m/z values at which to bisect spectrum
    :param out_csv_path: (string) Save path for csv of raw ion data
    :param img_out_path: (string) save path for ion images
    :param pv: Shared-state variable ("progress var") for usage with progress bar.
    :return: (list, obj) [(ndarray) average maldi image, (list, float) ions parsed]
    """
    # pv is 10 for imzmlparser, 100 for getionimage = 110
    logging.info('parsing ImzML file')
    logging.info('ions: '+ str(ions))
    p = ImzMLParser(MALDI_path)
    pv.set(pv.get()+10)
    xmax, xmin, ymax, ymin = getmaxmin(p)

    tolerances = [0.3]*len(ions)

    out_csv_path = project_dir + '/raw ion data.csv'

    output_folder = project_dir + '/maldi images/'
    os.mkdir(output_folder)
    logging.info('reading spectra')
    t = time.time()
    avg = record_reader([xmin,ymin,ymax,xmax], p, output_folder, ions, tolerances, out_csv_path, pv)
    elapsed = time.time()-t
    logging.info('reading spectra elapsed time: '+str(elapsed)+' s')

    return [avg, ions]




