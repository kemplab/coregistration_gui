import tkinter as tk
from tkinter.ttk import Progressbar
from tkinter.font import Font, BOLD

import Page
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

"""
@author Tanya Roysam
roysamt@gatech.edu
Output_Page.py
Child class of Page (Page.py). Specialized Page type for displaying output of pipeline steps or intermediate pages.
Not for input pages.
"""


class Output_Page(Page.Page):
    """
    Inherited from Page.
    Specialized Page type for displaying output of pipeline steps or intermediate pages.
    Not for input pages.

    * Individual "--gen" methods called ONLY during initial construction of page
    """
    FROM_COREG = 0
    FROM_IMZML = 1
    FROM_SEG = 2
    FROM_MET = 3
    def __init__(self, root, title, funcnum, ls, method=None):
        # parent constructor (from page class)-creating a "blank" page with no inputs and title title
        super().__init__(root, title)

        self.func_num = funcnum

        # creating new itemframe for objects
        self.item_frame = tk.Frame(self)
        # attribute array of images for PIL
        self.images = []
        # Progress bar (for loading pages)
        self.progress = None


        # Depending on the part of the pipeline, a different page generator will be called
        if funcnum == 1:
            self.step_1(*ls)
        elif funcnum == 2:
            self.step_2(*ls)
        elif funcnum == 3:
            self.step_3(*ls)
        elif funcnum == 4 or funcnum == 7:
            self.step_4(ls[0], ls[1], method, ls[2], ls[3])
        elif funcnum == 5:
            self.step_5(*ls)
        elif funcnum == 6:
            self.step_6(*ls)
        elif funcnum == 10:
            self.start_page(method)
        elif funcnum == -1:
            self.error_page(ls[0])
        else:
            self.create_loading(*ls)

    def create_plt(self, images, labels, plot_points=False, img2=None, lbl2=None, X=None, Y=None, coord_plot=False):
        num = len(images)
        if plot_points:
            num += 1
        fig = Figure(figsize=(6,5))
        if coord_plot:
            axes = fig.subplots(nrows=1, ncols=num, sharex=True, sharey= True)
        else:
            axes = fig.subplots(nrows=1, ncols=num)
        if (num == 1):
            axes.imshow(images[0])
            axes.set_title(labels[0])
        else:
            for ind, img in enumerate(images):
                axes[ind].imshow(img)
                axes[ind].set_title(labels[ind])

        if plot_points:
            axes[num-1].imshow(img2)
            axes[num-1].set_title(lbl2)
            axes[num-1].scatter(Y, X, c='r',s=3)


        canvas = FigureCanvasTkAgg(fig, master=self.item_frame)
        canvas.draw()

        toolbar = NavigationToolbar2Tk(canvas, self.item_frame, pack_toolbar=False)
        toolbar.update()

        toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        widget = canvas.get_tk_widget()
        widget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    def start_page(self, method):

        tk.Label(self.item_frame, text='').pack() # TODO add some text here idk lmao

        font = Font(size=14, weight=BOLD)

        # these are local
        def cbcmd():
            method(Output_Page.FROM_COREG)
        def sbcmd():
            method(Output_Page.FROM_SEG)
        def ibcmd():
            method(Output_Page.FROM_IMZML)
        def mbcmd():
            method(Output_Page.FROM_MET)

        # TODO add segmentation?

        self.cb = tk.Button(self.item_frame, text='Coregistration', pady=25, width=25, font=font, command=cbcmd)
        self.sb = tk.Button(self.item_frame, text='Segmentation', pady=25, width=25, font=font, command=sbcmd)
        self.ib = tk.Button(self.item_frame, text='Parse ImzML', pady=25, width=25, font=font, command=ibcmd)
        self.mb = None#tk.Button(self.itemframe, text='Additional Metrics', pady=25, width=25, font=font, command=mbcmd)

        self.cb.pack(padx = 25)
        self.sb.pack(padx = 25)
        self.ib.pack(padx = 25)
        #self.mb.pack(padx=25)

        self.item_frame.pack()

        self.update_idletasks()



    def step_1(self, confimg, maldimg):
        """
        Output page for preprocessed images

        :param confimg: Initial preprocessed confocal image
        :param maldimg: Initial preprocessed MALDI MSI image
        :return: None
        """
        """
        # Labels for each image
        tk.Label(self.itemframe, text='Preprocessed Confocal').grid(row=0,column=0,padx=10,pady=10)
        tk.Label(self.itemframe, text='Preprocessed MALDI').grid(row=0,column=1,padx=10,pady=10)

        # Create PIL images
        ci = Image.fromarray((cv2.normalize(confimg,None,alpha=0,beta=255,norm_type=cv2.NORM_MINMAX)).astype(np.uint8))
        mi = Image.fromarray((cv2.normalize(maldimg,None,alpha=0,beta=255,norm_type=cv2.NORM_MINMAX)).astype(np.uint8))

        # Create PIL ImageTk from ci, mi
        confpil = ImageTk.PhotoImage(image=self.resize(ci,250))
        maldpil = ImageTk.PhotoImage(image=self.resize(mi,250))

        # Save in attribute so images will render
        self.images = [confpil,maldpil]
        """
        self.create_plt([confimg, maldimg], ['Preprocessed Confocal', 'Preprocessed MALDI'])

        # Display image
        """tk.Label(self.itemframe, image=self.images[0]).grid(row=1,column=0, padx=20,pady=20)
        tk.Label(self.itemframe, image=self.images[1]).grid(row=1,column=1, padx=20,pady=20)"""
        self.item_frame.pack()

        # update
        self.update_idletasks()

    def step_2(self, img, alignment, confocal3, maldi3, h, w, xopt):
        """
        Output page for initial alignment

        :param img: (ndarray, float) added image of initial alignment
        :param alignment: (float) initial alignment score
        :return: None
        """

        """# PIL ImageTk for display
        imgpil = ImageTk.PhotoImage(Image.fromarray((cv2.normalize(img,None,alpha=0,beta=255,norm_type=cv2.NORM_MINMAX)).astype(np.uint8)))

        # Saving in attribute so image will render
        self.images = [imgpil]

        # Display image and initial alignment label
        tk.Label(self.itemframe, image=self.images[0]).pack(padx=10, pady=10)"""

        self.create_plt([img], [''])
        tk.Label(self.item_frame, text=('Initial Alignment Score: ' + str(alignment))).pack(padx=10, pady=10)
        self.item_frame.pack()

        # update
        self.update_idletasks()

    def step_3(self, img, alignment, outstr, savedirs, cropconf, maldi, crop_params):
        """
        Output page for final alignment

        :param img: (ndarray, float) added image of final alignment
        :param alignment: (float) final alignment score
        :param outstr: (string) output string describing alignment params
        :param savedirs: (list, string) output directories for MALDI and confocal
        :return:
        """

        """# PIL ImageTK for display
        imgpil = ImageTk.PhotoImage(Image.fromarray((cv2.normalize(img,None,alpha=0,beta=255,norm_type=cv2.NORM_MINMAX)).astype(np.uint8)))

        # Saving in attribute so image renders
        self.images = [imgpil]

        # Displaying images and labels for scores, outstr, save dirs
        tk.Label(self.itemframe, image=self.images[0]).pack(padx=10, pady=10)"""
        self.create_plt([img], [''])

        tk.Label(self.item_frame, text=('Final Alignment Score: ' + str(alignment))).pack(padx=10, pady=10)
        tk.Label(self.item_frame, text=outstr).pack(padx=10, pady=10)
        tk.Label(self.item_frame, text=('Images saved in ' + savedirs[0] + ',\n ' + savedirs[1])).pack(padx=10, pady=10)
        self.item_frame.pack()

        # update
        self.update_idletasks()

    def step_4(self, img, nlabels, method, crconf, seg_save):
        """# PIL ImageTK for display
        rawimg = Image.fromarray((cv2.normalize(img,None,alpha=0,beta=255,norm_type=cv2.NORM_MINMAX)).astype(np.uint8))
        rawimg = self.resize(rawimg, 400)
        imgpil = ImageTk.PhotoImage(rawimg)

        # Saving in attribute so image renders
        self.images = [imgpil]

        # Displaying image
        tk.Label(self.itemframe, image=self.images[0]).pack(padx=10, pady=10)"""

        X = np.array(nlabels['centroid-0'])
        Y = np.array(nlabels['centroid-1'])

        med = np.median(crconf)
        std = np.std(crconf)

        # Thresholding outliers out
        # High-intensity outliers distort results so I truncate them
        thr = int(round(med + 1 * std))

        crconf[crconf>thr] = thr

        self.create_plt([img], [''], plot_points=True, img2=crconf, lbl2='Seeds on Original', X=X, Y=Y, coord_plot=True)

        tk.Label(self.item_frame, text='Output saved at {}'.format(seg_save)).pack(padx=10, pady=10)
        tk.Label(self.item_frame, text='To rerun segmentation, return to previous page, change parameters and hit Continue.\n'
                                      'To calculate additional cell network metrics, toggle the checkbox and hit Continue.').pack(padx=10, pady=10)

        tk.Checkbutton(self, text='Calculate Cell Network Metrics', command=method).pack(padx=10, pady=10)

        self.item_frame.pack()

        self.update_idletasks()

    def step_5(self, avg, ions):

        self.create_plt([avg], [''])
        tk.Label(self.item_frame, text='imzML parsing complete.' +
                                      '\nPlease exit and relaunch to rerun coregistration or parsing.').pack(padx=10,
                                                                                                             pady=10)

        self.item_frame.pack()

        self.update_idletasks()

    def step_6(self, cnm_out):
        tk.Label(self.item_frame, text='Metrics and visualizations saved at:\n' + str(cnm_out) +
                                '\nPlease exit and relaunch to rerun coregistration or parsing.').pack(padx=10, pady=10)
        self.item_frame.pack()
        self.update_idletasks()


    def create_loading(self, steps, intvar, det=True):
        """
        Generate generic loading page with progress bar

        :param steps: (int) Number of steps to completion
        :param intvar: Shared state variable for current steps completed
        :return: None
        """
        if det:
            m = 'determinate'
            self.progress = Progressbar(self.item_frame, orient=tk.HORIZONTAL, length=200, mode=m, maximum=steps,
                                        variable=intvar)
        else:
            m = 'indeterminate'
            self.progress = Progressbar(self.item_frame, orient=tk.HORIZONTAL, length=200, mode=m)
            self.progress.start()

        self.progress.pack()
        self.item_frame.pack()


    def error_page(self, ex):
        tk.Label(self.item_frame, text='An error occurred and the app shut down.\n'
                                      'Please check that your inputs are valid and relaunch.\n'
                                      'For more information, check the session log at <time of launch>.log\n'
                                      'in the folder you saved MALDI output.').pack(padx=10,pady=10)

        errstring = ''
        for i in ex.args:
            errstring += str(i)
            errstring += '\n'
        tk.Label(self.item_frame, text='Error type: ' + errstring).pack(padx=10, pady=10)
        self.item_frame.pack()


    def resize(self, img, maxs):
        """
        Resizing to standard (specified) maximum dimension for display

        :param img: (ndarray, float) image to be resized
        :param maxs: (int) maximum size of either dimension. Will resize along largest dimension to this maximum
        :return: im
        """
        # current size
        rw, rh = img.size

        # finding largest dim, finding scalar multiplier, then new width/height
        if rw > rh:
            scalar = rw / maxs
            nw = maxs
            nh = int(rh / scalar)
        else:
            scalar = rh / maxs
            nh = maxs
            nw = int(rw / scalar)

        return img.resize((nw, nh))
