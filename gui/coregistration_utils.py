# import libraries

import numpy as np
import cv2

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure

from scipy.optimize import dual_annealing
from scipy.spatial.distance import cdist


from skimage.filters import median, gaussian
from skimage.morphology import disk, label
from skimage.measure import regionprops_table
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from sklearn.cluster import KMeans, DBSCAN

from pyimzml.ImzMLParser import ImzMLParser,_bisect_spectrum

import pandas as pd
import os
import random
import math
import time

import logging
import errno

"""
@author Tanya Roysam 
roysamt@gatech.edu
coregistration_utils.py
Breaks up coregistration pipeline for usage in multiprocessing/threading.
First section is from coregistration_metrics.py
"""


def transform_image(X, theta, shiftX, shiftY, scale):
    (height, width) = X.shape[:2]
    d = (height * height * 0.25 + width * width * 0.25) ** 0.5
    beta = np.arctan(height / width)
    pad_size = np.abs(int(d * np.cos(beta - theta) - width / 2))
    vertical_pad = np.zeros((height, pad_size))
    horizontal_pad = np.zeros((pad_size, 2 * pad_size + width))
    X_pad = np.concatenate((vertical_pad, X, vertical_pad), axis=1)
    X_pad2 = np.concatenate((horizontal_pad, X_pad, horizontal_pad), axis=0)
    (height2, width2) = X_pad2.shape[:2]
    # rotation
    M = cv2.getRotationMatrix2D((width2 / 2, height2 / 2), np.degrees(theta), 1)
    Y = cv2.warpAffine(X_pad2, M, (width2, height2))
    # scaling
    Y = cv2.resize(Y, dsize=(int(width2 * scale), int(height2 * scale)))
    (height3, width3) = Y.shape[:2]
    padY = int((height3 - height) / 2)
    padX = int((width3 - width) / 2)
    # shift
    if height3 < height and width3 < width:
        output = np.zeros((height, width))
        output[-padY:-padY + height3, -padX:-padX + width3] = Y
        M = np.float32([[1, 0, shiftX], [0, 1, shiftY]])
        res = cv2.warpAffine(output, M, (width, height))
        return res
    elif height3 > height and width3 > width:
        M = np.float32([[1, 0, shiftX], [0, 1, shiftY]])
        res = cv2.warpAffine(Y, M, (width3, height3))
        return res[padY:padY + height, padX:padX + width]
    else:
        M = np.float32([[1, 0, shiftX], [0, 1, shiftY]])
        res = cv2.warpAffine(Y, M, (width3, height3))
        return cv2.resize(res, dsize=(width, height))


def mutual_information(p,X,Y,pv=None):
    global GLOBAL_MIN
    (height, width) = X.shape[:2]
    if int(p[3]*height)*int(p[3]*width) == 0: return 0
    transformed_Y = transform_image(Y,p[0],p[1],p[2],p[3])
    [top_maldi,bottom_maldi,left_maldi,right_maldi] = crop(transformed_Y)
    [top_confocal,bottom_confocal,left_confocal,right_confocal] = crop(X)
    top,bottom = max(top_maldi,top_confocal),min(bottom_maldi,bottom_confocal)
    left,right = max(left_maldi,left_confocal),min(right_maldi,right_confocal)
    transformed_Y = transformed_Y[top:bottom,left:right]
    X = X[top:bottom,left:right]
    # add a fine factor to account for the parts of the image that were placed outside of the initial field of view
    fine = np.sum(transformed_Y > 0)/np.sum(Y > 0)
    X = X.ravel()
    transformed_Y = transformed_Y.ravel()
    # calculating mutual information
    hgram, _, _ = np.histogram2d(X,transformed_Y,bins=20)
    pxy = hgram / float(np.sum(hgram))
    py = np.sum(pxy, axis=1)
    px = np.sum(pxy, axis=0)
    px_py = py[:, None] * px[None, :]
    nzs = pxy > 0


    if pv:
        pv.set(pv.get()+1)

    return -fine*fine*np.sum(pxy[nzs] * np.log(pxy[nzs] / px_py[nzs]))



def crop(image):
    # crop out stripes of zeros on the sides of the image
    shape = image.shape
    top, bottom, left, right = 0, shape[0], 0, shape[1]
    for i in range(shape[0]):
        s = np.sum(image[i, :])
        if s > 0:
            bottom = i
            if top == 0: top = i
    for i in range(shape[1]):
        s = np.sum(image[:, i])
        if s > 0:
            right = i
            if left == 0: left = i
    return top, bottom, left, right

def cleanup(cellprops):
    areas = cellprops['area']
    area = np.array(areas)
    area_std, area_mean = np.std(area), np.mean(area)

    # exclude unreasonably small and large cells
    area_mask = (area < (area_mean + 3 * area_std)) & (area > max(5, area_mean - 3 * area_std))
    cells = cellprops[area_mask]

    X = np.array(cells['centroid-0'])
    Y = np.array(cells['centroid-1'])
    perimeter = np.array(cells['perimeter'])
    eccentricity = np.array(cells['eccentricity'])
    return cells, np.array(X), np.array(Y), area, perimeter, eccentricity

def relabel_random(labels):
    # no more for loops :)
    m = np.max(labels)
    random_labels = np.linspace(1,m,m)
    random.shuffle(random_labels)
    seeds = np.where(labels!=0)
    # assign each seed to a random label
    labels[seeds] = random_labels
    return labels

def plot(image, name):
    image_std, image_mean = np.std(image), np.mean(image)
    image[np.where(image > (image_mean + 2 * image_std))] = image_mean + 2 * image_std
    colors_list = [(1, 1, 1), (0, 0, 1), (0, 1, 0), (1, 1, 0), (1, 0, 0)]
    cm = LinearSegmentedColormap.from_list('cmap_name', colors_list, N=100)
    fig = Figure()
    axes = fig.subplots()
    mp = axes.imshow(image, cmap=cm)
    fig.colorbar(mp, ax=axes)
    fig.savefig(name, dpi=1200, bbox_inches="tight")
    return fig


def overlay_cells(image, cells, name, dapi=None):
    h, w = image.shape
    output = np.zeros((h, w))
    n = len(cells)
    intensities = np.zeros(n)
    area = np.array(cells['area'])

    # iterate through each cell and find sum of pixel intensities
    for i in range(n):
        coords = np.array(cells.iloc[i].coords)
        intensities[i] = np.sum(image[coords[:, 0], coords[:, 1]])
    intensities /= area

    # set each coordinate corresponding to a cell to mean intensity value
    for i in range(n):
        coords = np.array(cells.iloc[i].coords)
        output[coords[:, 0], coords[:, 1]] = intensities[i]
    if dapi is not None:
        output[np.where(dapi != 0)] = output[np.where(dapi != 0)] / dapi[np.where(dapi != 0)]
    diagram = plot(output, name)
    return intensities, output, diagram


# visualize metrics, saving in the output folder provided through the "name" variable
def visualize(metric, cells, h, w, name):
    output = np.zeros((h, w))
    n = len(cells)
    for i in range(n):
        coords = cells.iloc[i].coords
        output[coords[:,0], coords[:,1]] = metric[i]
    plot(output, name)



# get average distance between cells' centers from a random rectangular sample in the image
def get_connection_length(X, Y):
    # select random region
    left = random.randrange(int(np.min(X)), int(0.75 * np.max(X)))
    top = random.randrange(int(np.min(Y)), int(0.75 * np.max(Y)))
    rand_ind = np.where(np.logical_and(np.logical_and(X > left, X < (left + 0.25 * np.max(X))),
                                       np.logical_and(Y > top, Y < (top + 0.25 * np.max(Y)))))
    X = X[rand_ind]
    Y = Y[rand_ind]
    n = len(X)
    connections = np.zeros(n)
    for i in range(n):
        d = np.max(X)
        for j in range(n):
            if i == j: continue
            dist = ((X[j] - X[i]) ** 2 + (Y[j] - Y[i]) ** 2) ** 0.5
            if dist < d: d = dist
        connections[i] = d
    mu = np.mean(connections)
    sigma = 3 * np.std(connections)
    main_body = np.logical_and(connections < (mu + sigma), connections > (mu - sigma))
    connections = connections[main_body]
    l = np.max(connections)
    return l


def edge_indiv(selfX, selfY, neighX, neighY):
    # find components of distance vector between center cell and its neighbors: nx1 arrays
    dist_vec_y = selfY-neighY
    dist_vec_x = selfX-neighX

    # find angles of distance vectors
    angles = np.arctan2(dist_vec_y,dist_vec_x)
    angles[angles<0] = np.pi * 2 + angles[angles<0]

    # sort angles
    sorted_angles = np.sort(angles) # default quicksort
    if len(sorted_angles)==1: return True
    difference = np.max(np.abs(np.diff(sorted_angles))) # Does not include beginning and last
    # Include distance between smallest and largest angle
    last_diff = np.pi * 2 - (sorted_angles[len(sorted_angles)-1] - sorted_angles[0])

    # true or false
    return np.max([difference, last_diff]) >= (np.pi / 2)

def edge(X, Y, global_neigh_region, has_connections, edge_margin):
    # only cells that have connections
    nz_inds = np.where(has_connections)[0]

    # all X/Y indices of neighboring cells for each cell. list of arrays
    global_neigh_region_x = [X[locs] for locs in global_neigh_region]
    global_neigh_region_y = [Y[locs] for locs in global_neigh_region]

    edge1 = np.ones(len(X))

    nzX = X[nz_inds]
    nzY = Y[nz_inds]

    gnrx_nz = [global_neigh_region_x[i] for i in nz_inds]
    gnry_nz = [global_neigh_region_y[i] for i in nz_inds]

    # for cells that have neighbors: edge cell calculation
    edge_nz_1 = np.array([edge_indiv(nzX[i],nzY[i],gnrx_nz[i],gnry_nz[i]) for i in range(len(nz_inds))])

    edge1[has_connections] = edge_nz_1

    # cells that are edge_margin from the edge of the image should not be considered edge cells
    bordermask = ((X < edge_margin) | (Y < edge_margin) | (X > max(X) - edge_margin) | (
                Y > max(Y) - edge_margin)) & has_connections
    edge1[bordermask] = 0

    return edge1.astype(int)


# calculating distance from the edge of the colony for a given cell
def edgeDistance(X, Y, edge, i):
    distances = ((X[np.where(edge == 1)] - X[i]) ** 2 + (Y[np.where(edge == 1)] - Y[i]) ** 2) ** 0.5
    return np.min(distances)


def ClusterImage(Z, k, X, Y, epsilon, minclust):
    if np.sum(Z) == 0: return np.zeros(len(Z)), np.ones(len(Z))
    k = int(round(k))
    kmeans = KMeans(n_clusters=k, random_state=0).fit(Z.reshape(-1, 1))
    idx = kmeans.labels_
    l = len(Z)
    cluster_sizes = np.zeros(l)
    cluster_avg = np.zeros(l)
    nums = np.arange(0, l)
    randomness = 0
    for i in range(k):
        indices = nums[np.where(idx == i)]
        x = X[indices]
        y = Y[indices]
        z = Z[indices]
        coeffs = nums[indices]
        clustering = DBSCAN(eps=epsilon, min_samples=minclust).fit(np.array([x, y]).reshape(-1, 2))
        idx2 = clustering.labels_
        dump = 0
        sizes = []
        avgs = []
        for j in range(len(idx2)):
            if idx2[j] == -1:
                sizes.append(1)
                if z[j] == 0:
                    avgs.append(0.01)
                else:
                    avgs.append(z[j])
                dump += 1
            else:
                sizes.append(np.sum(idx2 == idx2[j]))
                if np.sum(z[np.where(idx2 == idx2[j])]) == 0:
                    avgs.append(0.01)
                else:
                    avgs.append(np.sum(z[np.where(idx2 == idx2[j])]) / np.sum(idx2 == idx2[j]))
        cluster_sizes[coeffs] = sizes
        cluster_avg[coeffs] = avgs
        randomness = randomness + np.max(idx2) + 2 * dump
        nums[np.where(idx == i)] = idx2 + np.max(nums) * np.ones(len(idx2)) - np.array(idx2 == 0, dtype=int) * np.max(
            nums)

    randomness = randomness / len(Z)

    return cluster_sizes, cluster_avg


def getmaxmin(p):
    xmax, xmin, ymax, ymin = 0, 600000, 0, 600000
    for i, (x, y, z) in enumerate(p.coordinates):
        if x > xmax: xmax = x
        if x < xmin: xmin = x
        if y > ymax: ymax = y
        if y < ymin: ymin = y
    return xmax, xmin, ymax, ymin


def getionimage(borders,p, mz_values, tolerances, z=1, reduce_func=sum,dim=1,pv=None):
    coordsX = []
    coordsY = []
    im = np.zeros((dim, borders[2] - borders[1] + 1, borders[3] - borders[0] + 1))
    sum_im = np.zeros((borders[2] - borders[1] + 1, borders[3] - borders[0] + 1))

    total = len(p.coordinates)
    percent = int(total/100)


    for i, (x, y, z_) in enumerate(p.coordinates):
        if z_ == z and x >= borders[0] and y >= borders[1] and x <= borders[3] and y <= borders[2]:
            mzs, ints = p.getspectrum(i)
            sum_im[y - borders[1] - 1, x - borders[0] - 1] = np.sum(ints)
            coordsX.append(x - borders[0])
            coordsY.append(y - borders[1])

            for index, mz_value in enumerate(mz_values):
                min_i, max_i = _bisect_spectrum(mzs, mz_value, tolerances[index])
                integral_signal = reduce_func(ints[min_i:max_i + 1])
                im[index, y - borders[1] - 1, x - borders[0] - 1] = integral_signal
        if i%percent==0: pv.set(pv.get()+1)

    return im, sum_im, coordsX, coordsY

def record_reader(borders,p,MALDI_output,mz_values,tolerances, csv_output, pv):

    img, sum_im, coordsX, coordsY = getionimage(borders, p, mz_values, tolerances, z=1, reduce_func=sum,
                                                dim=len(mz_values), pv=pv)

    average = np.sum(img, axis=0) / len(mz_values)
    plt.imsave(os.path.join(MALDI_output, "average.png"), average)

    for index, values in enumerate(img):
        a = values
        exten = str(mz_values[index])
        exten.replace(".", "_")
        plt.imsave("{}//MALDI__{}.png".format(MALDI_output, exten), a)

    headings = {'X': coordsX, 'Y': coordsY}
    for index, values in enumerate(img):
        vals = 1000 * values[coordsY, coordsX]  / sum_im[coordsY, coordsX] #sum of ions change
        headings.update({str(mz_values[index]): vals})

    for index, values in enumerate(img):
        vals = 1000 * values[coordsY, coordsX]  # / sum_im[coordsY, coordsX]
        headings.update({'raw (non-normalized)'+str(mz_values[index]): vals})

    sum_vals = sum_im[coordsY, coordsX]
    headings.update({'spectrum integral': sum_vals})

    df = pd.DataFrame(headings)
    df.to_csv(csv_output)

    return average


def center(img, inh, inw):
    h,w = img.shape

    #h shift
    th = int(math.floor((inh-h)/2))

    #w shift
    lw = int(math.floor((inw-w)/2))

    centered = np.zeros((inh, inw))
    centered[th:th+h, lw:lw+w] = img

    return centered

def ASF(image, filter_size, normalize=False):
    """
    Alternating Sequential Filter

    :param image: input grayscale image (width,height)
    :param filter_size: list of filter size
    :param normalize: normalize the image to dtype
    """
    background = np.copy(image)
    for sz in filter_size:
        # apply morphological opening with the defined structuring element
        selem = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (sz, sz))
        background = cv2.morphologyEx(background, cv2.MORPH_OPEN, selem)
        background = cv2.morphologyEx(background, cv2.MORPH_CLOSE, selem)

    # subtract the background from image
    image = cv2.subtract(image, background)

    # normalize image to dtype
    if normalize:
        return cv2.normalize(image, None, norm_type=cv2.NORM_MINMAX, dtype=-1)
    else:
        return image

def getimagedimsfromcsv(file):
    df = pd.read_csv(file)
    maxes = df.max(axis=0)
    ydim = maxes['Y'].astype(int)
    xdim = maxes['X'].astype(int)
    return ydim,xdim

def readimagefromcsv(mz_value,file,ydim,xdim):
    df = pd.read_csv(file)

    x = df['X'].to_numpy()
    y = df['Y'].to_numpy()
    val = df[str(mz_value)].to_numpy()

    img = np.zeros((int(ydim)+1,int(xdim)+1))

    try:
        img[y,x] = val
    except IndexError:
        raise IndexError('Check that every row in raw MALDI data csv has associated non-negative integer X and Y indices.')

    return img

def readimagefromdf(mz_value, df, ydim, xdim):
    x = df['X'].to_numpy()
    y = df['Y'].to_numpy()
    val = df[str(mz_value)].to_numpy()

    img = np.zeros((int(ydim) + 1, int(xdim) + 1))

    try:
        img[y,x] = val
    except IndexError:
        raise IndexError('Check that every row in raw MALDI data csv has associated non-negative integer X and Y indices.')

    return img
