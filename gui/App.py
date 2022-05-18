import Page as pg
import Output_Page as pf
import coregistration_pipeline
import os
import Navigator as pv
import re
import time
import logging
from getpass import getuser
from stat import S_IROTH as read_only
from sys import version, platform

"""
@author Tanya Roysam
roysamt@gatech.edu
App.py
Handles function of overall GUI. 
"""


class App(pv.Navigator):
    """
    Child class of Navigator. Overall navigation object that contains Page objects.
    Handles function of overall GUI.
    * Navigation from page to page
    * Checking and fetching user input
    * Initiating secondary threads and handling output of threads
    """
    # different pathways
    FROM_COREG = 0
    FROM_IMZML = 1
    FROM_SEG = 2
    FROM_COREG2 = 3
    FROM_SEG2 = 4

    # steps
    START = 0
    COREG_IP = 1
    PREPROC = 2
    INIT_ALIGN = 3
    DUAL_ANNEAL = 4
    SEG_IN = 5
    SEG = 7
    IMZML_IN = 8
    IMZML = 9
    REPEAT_SEG = 10
    CELL_METRICS_IP = 11
    CELL_METRICS = 12
    SEG_IP = 13
    REPEAT_SEG_IP = 14
    FINISH = 15

    # content of each pathway
    # TODO: swap each coregistration/segmentation order so that repeating segmentation is not default?
    FC_ORDER = [COREG_IP, PREPROC, INIT_ALIGN, DUAL_ANNEAL, SEG_IN, SEG, REPEAT_SEG]
    FC_ORDER_2 = [COREG_IP, PREPROC, INIT_ALIGN, DUAL_ANNEAL, SEG_IN, SEG, CELL_METRICS_IP, CELL_METRICS, FINISH]
    FI_ORDER = [IMZML_IN, IMZML, FINISH]
    FS_ORDER = [SEG_IP, SEG, REPEAT_SEG_IP]
    FS_ORDER_2 = [SEG_IP, SEG, CELL_METRICS_IP, CELL_METRICS, FINISH]

    def __init__(self, root):
        """
        Constructor for overall app

        :param root: (tk) Master window
        """


        methods = [coregistration_pipeline.preprocessing, coregistration_pipeline.initial_alignment,
                   coregistration_pipeline.final_alignment, coregistration_pipeline.segment,
                   coregistration_pipeline.imzml_parse, coregistration_pipeline.cell_network,
                   coregistration_pipeline.pre_segmentation]

        output_titles = ['Preprocessed', 'Initial Alignment', 'Final Alignment', 'Segmented Confocal',
                         'Average MALDI image of selected ions', 'Cell network metrics', 'Segmented Confocal']

        arg_names = [['MALDI MSI average image', 'Confocal image', 'pv', 'minsz',
                     'maxsz', 'Remove background? Only recommended for sparsely scattered cells, not dense colonies.'],

                    ['confocal2','maldi2','xopt','pv','Scale of the MALDI image (um/px)',
                     'Scale of the confocal image (um/px)'],

                    ['confocal3', 'maldi3', 'h', 'w', 'Project folder',
                     'MALDI MSI average image', 'Confocal image', 'xopt', 'Maximum global annealment iterations',
                     'Range of angle', 'Range of x-shift', 'Range of y-shift', 'Range of scale', 'pv',
                     'MALDI raw ion intensity data csv', 'ions'],

                    ['cropped_confocal', 'Minimum radius of cell to detect (um).',
                     'Maximum radius of cell to detect (um).', 'Project folder',
                     'Colony is made up of bright, separated, scattered cells',
                     'Automatically remove seeds not within colony (recommended for clustered colonies only)',
                     'Number of thresholds. Increase if cells are not detected.',
                     'For images with high-intensity noise: Threshold (in standard deviations above mean).',
                     'resized_maldis', 'ions', 'crop_params', 'Confocal image', 'pv'],

                    ['ImzML file', 'ions', 'Project folder', 'pv'],

                    ['seg_save', 'Cell region radius (int, in cell-cell connections) e.g. 3',
                     'Intensity classes for k-means clustering',
                     'Minimum size of k-means cluster', 'connection',
                     'Project folder', 'edge margin', 'cropped_confocal', 'labels df',
                     'Visualize metrics? Will take time.', 'pv'],


                     ['MALDI MSI average image', 'Confocal image', 'xopt', 'MALDI raw ion intensity data csv',
                      'Minimum radius of cell to detect (um).', 'Maximum radius of cell to detect (um).',
                      'Project folder', 'Colony is made up of bright, separated, scattered cells',
                      'Automatically remove seeds not within colony (recommended for clustered colonies only)',
                      'Number of thresholds. Increase if cells are not detected.',
                      'For images with high-intensity noise: Threshold (in standard deviations above mean).',
                      'ions', 'Confocal image', 'pv']
                    ]

        # create start page
        self.start_page = pf.Output_Page(root, title='Start', funcnum=10, ls=None,
                                         method=self.start)
        page_list = [self.start_page]

        # parent constructor called
        super(App, self).__init__(root, page_list, methods, output_titles, arg_names)

        self.cont_button.config(state='disabled', command=self.next_step)

        self.methods = methods
        self.current_step = App.START

        self.pathway = None

        self.output_titles = output_titles

    def start(self, pathway):
        """
        Select pathway on start screen

        :param pathway: Pathway number class var
        :return: None
        """

        self.pathway = pathway
        self.cont_button.config(state='normal')

        if pathway==App.FROM_COREG:
            self.start_coregistration()
        elif pathway==App.FROM_SEG:
            self.start_segmentation()
        else:
            self.start_imzml()

    def parse_ions(self, fn):
        """
        Parse .txt of ions to list of floats
        :param fn: Name of .txt
        :return: (list) ions
        """
        try:
            fin = open(fn, 'r')
        except FileNotFoundError as e:
            self.error_handling(e)
            return -1
        strs = fin.readlines()
        fin.close()
        ion_str = ''
        for s in strs:
            ion_str += s
            ion_str += ','
        ion_str = re.split('[^0-9^.]', ion_str)

        ions = [float(x) for x in ion_str if x != '']

        return ions


    def special_ips_1(self):
        """
        Interpret certain special inputs, step 1 coregistration
        :return: None
        """

        confocal_pixel_size_um = self.inputs['Scale of the confocal image (um/px)']
        ions = self.parse_ions(self.inputs['Comma-separated list of ions (e.g. 883.576,885.573,887.57)'])


        init_rot = self.inputs['Initial MALDI rotation']
        init_xshift = int(self.inputs['Initial MALDI x-shift'])
        init_yshift = int(self.inputs['Initial MALDI y-shift'])
        init_scale = self.inputs['Initial MALDI scale']

        d = {'minsz':int(round(self.inputs['Minimum cell width (um)']/confocal_pixel_size_um)),
             'maxsz':int(round(self.inputs['Maximum cell width (um)']/confocal_pixel_size_um)),
             'ions': ions,
             'xopt':[init_rot, init_xshift, init_yshift, init_scale],
             'Maximum global annealment iterations':int(self.inputs['Maximum global annealment iterations']),
             'Range of x-shift':int(self.inputs['Range of x-shift']),
             'Range of y-shift':int(self.inputs['Range of y-shift'])
             }
        t = time.asctime()
        fn = self.inputs['Project folder'] + '/' + t.replace(':', '_') + '.log'
        logging.basicConfig(filename=fn, level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
        try:
            os.chmod(fn, read_only)
        except:
            logging.warning('Failed to make logfile read-only. Check permissions/ownership of destination directory')
        try:
            logging.info('User: ' + getuser())
            logging.info('Python {} on {}'.format(version, platform))
        except:
            pass

        self.inputs.update(d)

    def special_ips_2(self):
        """
        Interpret certain special inputs, segmentation input
        :return: None
        """
        conf_res = self.inputs['Scale of the confocal image (um/px)']
        low_rad = int(round(self.inputs['Minimum radius of cell to detect (um).'] / conf_res))
        high_rad = int(round(self.inputs['Maximum radius of cell to detect (um).'] / conf_res))
        num_thresh = int(round(self.inputs['Number of thresholds. Increase if cells are not detected.']))
        d = {'Minimum radius of cell to detect (um).':low_rad,
             'Maximum radius of cell to detect (um).':high_rad,
             'Number of thresholds. Increase if cells are not detected.':num_thresh
             }
        self.inputs.update(d)

    def special_ips_3(self):
        """
        Interpret certain special inputs, ImzML parsing
        :return: None
        """
        ions = self.parse_ions(self.inputs['Comma-separated list of ions (e.g. 883.576,885.573,887.57)'])

        d = {'ions':ions}

        t = time.asctime()
        fn = self.inputs['Project folder'] + '/' + t.replace(':', '_') + '.log'
        logging.basicConfig(filename=fn, level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
        try:
            os.chmod(fn, read_only)
        except:
            logging.warning('Failed to make logfile read-only. Check permissions/ownership of destination directory')
        try:
            logging.info('User: ' + getuser())
            logging.info('Python {} on {}'.format(version, platform))
        except:
            pass

        self.inputs.update(d)

    def special_ips_4(self):
        """
        Interpret certain special inputs, cell network metrics input
        :return: None
        """
        conf_res = self.inputs['Scale of the confocal image (um/px)']
        connection = int(round(self.inputs['Length of cell-cell connection in um (-1 to generate automatically)']/conf_res))
        edge_margin = int(round(self.inputs['Outer margin of image to ignore when detecting edge cells, in pixels']))

        d = {'connection':connection, 'edge margin':edge_margin}
        self.inputs.update(d)

    def special_ips_5(self):
        """
        Interpret certain special inputs, initial segmentation input
        :return:
        """
        init_rot = self.inputs['MALDI rotation']
        init_xshift =self.inputs['MALDI x-shift']
        init_yshift = self.inputs['MALDI y-shift']
        init_scale = self.inputs['MALDI scale']
        conf_res = self.inputs['Scale of the confocal image (um/px)']

        low_rad = int(round(self.inputs['Minimum radius of cell to detect (um).'] / conf_res))
        high_rad = int(round(self.inputs['Maximum radius of cell to detect (um).'] / conf_res))
        num_thresh = int(round(self.inputs['Number of thresholds. Increase if cells are not detected.']))

        ions = self.parse_ions(self.inputs['Comma-separated list of ions (e.g. 883.576,885.573,887.57)'])

        #cropped_confocal, confocal_crop_params, resized_maldis = multipipeline.pre_seg(maldi_path, conf_path, xopt, ions, csv)

        #d = {'cropped_confocal':cropped_confocal, 'crop_params':confocal_crop_params, 'resized_maldis':resized_maldis,
        d = {'ions':ions, 'Minimum radius of cell to detect (um).':low_rad,
             'Maximum radius of cell to detect (um).':high_rad,
             'Number of thresholds. Increase if cells are not detected.':num_thresh, 'xopt':[init_rot, init_xshift, init_yshift, init_scale]}

        t = time.asctime()
        fn = self.inputs['Project folder'] + '/' + t.replace(':', '_') + '.log'
        logging.basicConfig(filename=fn, level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
        try:
            os.chmod(fn, read_only)
        except:
            logging.warning('Failed to make logfile read-only. Check permissions/ownership of destination directory')
        try:
            logging.info('User: ' + getuser())
            logging.info('Python {} on {}'.format(version, platform))
        except:
            pass

        self.inputs.update(d)

    def toggle_path(self):
        """
        For segmentation PF. Toggles path from repeating segmentation to moving on and vice versa
        :return: None
        """
        if self.current_step == App.SEG:
            if self.pathway == App.FROM_COREG2:
                self.pathway = App.FROM_COREG
            elif self.pathway == App.FROM_COREG:
                self.pathway = App.FROM_COREG2
            elif self.pathway == App.FROM_SEG2:
                self.pathway = App.FROM_SEG
            elif self.pathway == App.FROM_SEG:
                self.pathway = App.FROM_SEG2

    def next_step(self):
        """
        Bound to continue.
        Move on to next step based on pathway in class variable
        :return: None
        """

        # Finding next step
        try:
            if self.pathway==App.FROM_COREG:
                next_step = App.FC_ORDER[App.FC_ORDER.index(self.current_step)+1]
            elif self.pathway==App.FROM_SEG:
                next_step = App.FS_ORDER[App.FS_ORDER.index(self.current_step)+1]
            elif self.pathway==App.FROM_COREG2:
                next_step = App.FC_ORDER_2[App.FC_ORDER_2.index(self.current_step) + 1]
            elif self.pathway==App.FROM_SEG2:
                next_step = App.FS_ORDER_2[App.FS_ORDER_2.index(self.current_step) + 1]
            else:
                next_step = App.FI_ORDER[App.FI_ORDER.index(self.current_step)+1]
        except IndexError:
            print('index error')
            next_step = App.FINISH

        # Storing outputs in current_out and updating inputs dict
        if self.current_step == App.PREPROC:
            confocal2, maldi2 = self.current_out
            self.inputs.update({'confocal2': confocal2, 'maldi2': maldi2})
        elif self.current_step == App.INIT_ALIGN:
            added_image, sc, confocal3, maldi3, h, w, xopt = self.current_out
            self.inputs.update({'confocal3': confocal3, 'maldi3': maldi3, 'h': h, 'w': w, 'xopt': xopt})
        elif self.current_step == App.DUAL_ANNEAL:
            finimg, finalsc, outstr, outpaths, cropped_confocal, resized_maldis, confocal_crop_params = self.current_out
            self.inputs.update({'cropped_confocal':cropped_confocal, 'resized_maldis':resized_maldis,
                                'crop_params':confocal_crop_params})

        elif self.current_step == App.SEG and (self.pathway == App.FROM_SEG or self.pathway == App.FROM_SEG2):
            nuclei_labels, cells, cropped_confocal, seg_save, confocal_crop_params, resized_maldis = self.current_out
            self.inputs.update({'cropped_confocal': cropped_confocal, 'crop_params': confocal_crop_params,
                                'resized_maldis': resized_maldis})

        if next_step==App.REPEAT_SEG:
            print('repeating segmentation')
            next_step = App.SEG
            self.current_step = self.SEG_IN
            self.del_page(len(self.page_list))
        elif next_step==App.REPEAT_SEG_IP:
            print('delete page, next step is segmentation')
            next_step = App.SEG
            self.current_step = self.SEG_IP
            self.del_page(len(self.page_list))


        # Handling each step.
        if next_step == App.PREPROC:

            yes = self.check([1])#,2,3])
            if yes:
                self.cont_text.set(' ')
                self.special_ips_1()
                args = self.get_args(1)
                self.current_step = App.PREPROC
                self.create_thread(args, self.methods[0], 8, 'Preprocessing...',1)
            else:
                self.cont_text.set('Check that all inputs are complete and valid.')

        elif next_step == App.INIT_ALIGN:

            args = self.get_args(2)
            self.current_step = App.INIT_ALIGN
            self.create_thread(args, self.methods[1], 1, None,2)

        elif next_step == App.DUAL_ANNEAL:

            args=self.get_args(3)
            self.current_step = App.DUAL_ANNEAL
            iters = self.inputs['Maximum global annealment iterations']
            self.create_thread(args, self.methods[2], iters*10, 'Aligning..',3)

        elif next_step == App.SEG_IN:

            self.current_step = App.SEG_IN
            self.segmentation_input()



        elif next_step == App.SEG:
            if self.pathway == App.FROM_COREG or self.pathway == App.FROM_COREG2:
                yes = self.check([5])
            else:
                yes = self.check([1])
            if yes:
                if self.pathway == App.FROM_COREG or self.pathway == App.FROM_COREG2:
                    self.special_ips_2()

                    self.cont_text.set(' ')

                    args = self.get_args(4)
                    self.current_step = App.SEG
                    nt = self.inputs['Number of thresholds. Increase if cells are not detected.']
                    high = self.inputs['Maximum radius of cell to detect (um).']
                    low = self.inputs['Minimum radius of cell to detect (um).']
                    ions = self.inputs['ions']
                    if not self.inputs['Colony is made up of bright, separated, scattered cells']:
                        # if normal segmentation
                        pvnum = 14 + 2 * nt + high - low + len(ions)
                    else:
                        # scattered cells
                        pvnum = 14 + high - low + len(ions)

                    self.create_thread(args, self.methods[3], pvnum, 'Segmenting and Overlaying...', 4,
                                       method_cbox=self.toggle_path)

                else:
                    self.special_ips_5()
                    self.cont_text.set(' ')

                    args = self.get_args(7)
                    self.current_step = App.SEG

                    nt = self.inputs['Number of thresholds. Increase if cells are not detected.']
                    high = self.inputs['Maximum radius of cell to detect (um).']
                    low = self.inputs['Minimum radius of cell to detect (um).']
                    ions = self.inputs['ions']

                    if not self.inputs['Colony is made up of bright, separated, scattered cells']:
                        # if normal segmentation
                        pvnum = 14 + 2 * nt + high - low + len(ions)
                    else:
                        # scattered cells
                        pvnum = 14 + high - low + len(ions)

                    self.create_thread(args, self.methods[6], pvnum, 'Segmenting and Overlaying...', 7,
                                       method_cbox=self.toggle_path)



            else:
                self.cont_text.set('Check that all inputs are complete and valid.')

        elif next_step == App.CELL_METRICS_IP:
            img, labels_df, cc, s_s, _, _ = self.current_out
            self.inputs.update({'labels df': labels_df, 'seg_save':s_s})
            self.current_step = App.CELL_METRICS_IP
            self.metrics_input()

        elif next_step == App.CELL_METRICS:
            if self.pathway == App.FROM_COREG or self.pathway == App.FROM_COREG2:
                yes = self.check([7])
            else:
                yes = self.check([3])
            if yes:
                self.cont_text.set(' ')
                self.special_ips_4()
                args = self.get_args(6)
                self.current_step = App.CELL_METRICS
                self.create_thread(args, self.methods[5], 8, 'Calculating cell network metrics...', 6)
            else:
                self.cont_text.set('Check that all inputs are complete and valid.')


        elif next_step == App.IMZML:
            yes = self.check([1])
            if yes:
                self.cont_text.set(' ')
                self.special_ips_3()
                args = self.get_args(5)
                self.current_step = App.IMZML
                self.create_thread(args, self.methods[4], 110, 'Parsing...', 5)
            else:
                self.cont_text.set('Check that all inputs are complete and valid.')

        elif next_step == App.FINISH:
            self.cont_button.config(state='disabled')

        else:
            self.cont_button.config(state='disabled')


    def start_imzml(self):
        """
        Start parsing ImzML
        :return: None
        """
        imzml_ip = self.imzml_input()
        self.add_page(imzml_ip)
        self.del_page(1)
        self.pathway = App.FROM_IMZML
        self.cont_button.config(command=self.next_step, state='normal')
        self.current_step = App.IMZML_IN

    def start_coregistration(self):
        """
        Start coregistration beginning with input
        :return: None
        """
        img_paths = ['MALDI MSI average image', 'Confocal image']
        misc_paths = ['MALDI raw ion intensity data csv', 'Comma-separated list of ions (e.g. 883.576,885.573,887.57)']
        folders = ['Project folder']
        values = {'Scale of the MALDI image (um/px)': 'float',
                   'Scale of the confocal image (um/px)': 'float',
                   'Minimum cell width (um)': 'float', 'Maximum cell width (um)': 'float'}

        defaults = {'Initial MALDI rotation': 0, 'Initial MALDI x-shift': 0, 'Initial MALDI y-shift': 0,
                   'Initial MALDI scale': 1, 'Maximum global annealment iterations': 1000, 'Range of angle': .025,
                   'Range of x-shift': 100, 'Range of y-shift': 100, 'Range of scale': .025}
        checkboxes = ['Remove background? Only recommended for sparsely scattered cells, not dense colonies.']
        page1 = pg.Page(self, 'Input the following parameters, upload the specified files below,\n '
                              'and create a project directory to save output.\n'
                              'Please upload a .txt file containing the desired ions.\n'
                              'Ions must be in specified csv.\n', values, [], img_paths, 'Advanced Settings', defaults,
                        checkboxes=checkboxes, misc_files=misc_paths, folders=folders, scroll=True)
        #page2 = pg.Page(self,'Create and select a project folder.\nAll output from this session will be saved here.',
                        #{}, [], [], folders=folderlist)
        #page3 = pg.Page(self, 'Fill in the following values.', valdict, [], [], 'Advanced Settings', defdict,
                        #checkdict)

        self.current_step = App.COREG_IP

        self.add_page(page1)
        #self.addpage(page2)
        #self.addpage(page3)

        self.del_page(1)

        self.pathway = App.FROM_COREG

        self.cont_button.config(command=self.next_step, state='normal')


    def segmentation_input(self):
        """
        Start segmentation with inputs
        :return: None
        """
        val_dict = {'Minimum radius of cell to detect (um).': 'float',
                    'Maximum radius of cell to detect (um).': 'float'}
        check_list = ['Colony is made up of bright, separated, scattered cells',
                      'Automatically remove seeds not within colony (recommended for clustered colonies only)']
        advanced_dict = {'Number of thresholds. Increase if cells are not detected.': 10,
                         'For images with high-intensity noise: Threshold (in standard deviations above mean).': 1.25}

        seg_ip = pg.Page(self, "Select the following parameters for segmentation.\n"
                               "To detect more or fewer cells, change the thresholds in the advanced settings.",
                         field_types=val_dict, default_text='Advanced Settings', default_fields=advanced_dict,
                         checkboxes=check_list, scroll=True)

        self.add_page(seg_ip)

    def imzml_input(self):
        """
        Generate imzml input page
        :return: (Page) ImzML input page object
        """
        misc_paths = ['ImzML file', 'Comma-separated list of ions (e.g. 883.576,885.573,887.57)']
        folder = ['Project folder']
        imzml_ip = pg.Page(self, "Select the following input files and create a project directory to save output.\n"
                                 "The corresponding .ibd file MUST have the same filename as the specified ImzML file and MUST\n"
                                 "be in the same directory.",
                           misc_files=misc_paths, folders=folder, scroll=True)
        return imzml_ip


    def metrics_input(self):
        """
        Create and add metrics input page
        :return: None
        """
        val_dict = {'Cell region radius (int, in cell-cell connections) e.g. 3': 'float'}
        advanced_dict = {'Intensity classes for k-means clustering': 3,
                         'Minimum size of k-means cluster': 5,
                         'Length of cell-cell connection in um (-1 to generate automatically)':-1,
                         'Outer margin of image to ignore when detecting edge cells, in pixels':50}
        checks = ['Visualize metrics? Will take time.']
        title = 'Select the following inputs to calculate cell network metrics.\n' \
                'Note that the cell region radius\n' \
                'is measured in individual cell-cell connection lengths,\n' \
                'which can be changed in the settings.'
        metrics_ip = pg.Page(self, title, field_types=val_dict,  default_text='Settings', default_fields=advanced_dict,
                             checkboxes=checks, scroll=True)

        self.add_page(metrics_ip)

    def start_segmentation(self):
        img_paths = ['MALDI MSI average image', 'Confocal image']
        misc_paths = ['MALDI raw ion intensity data csv', 'Comma-separated list of ions (e.g. 883.576,885.573,887.57)']
        folder_paths = ['Project folder']

        val_dict = {'MALDI rotation': 'float', 'MALDI x-shift': 'float', 'MALDI y-shift': 'float',
                    'MALDI scale': 'float','Minimum radius of cell to detect (um).': 'float',
                    'Maximum radius of cell to detect (um).': 'float', 'Scale of the MALDI image (um/px)': 'float',
                   'Scale of the confocal image (um/px)': 'float'}

        check_list = ['Colony is made up of bright, separated, scattered cells',
                      'Automatically remove seeds not within colony (recommended for clustered colonies only)']
        advanced_dict = {'Number of thresholds. Increase if cells are not detected.': 10,
                         'For images with high-intensity noise: Threshold (in standard deviations above mean).': .75}

        seg_ip = pg.Page(self, "Select the following parameters for segmentation.\n"
                               "Make sure to input the MALDI alignment parameters\n"
                               "exactly as they were output in the prior coregistration run.\n"
                               "To detect more or fewer cells, change the thresholds in the advanced settings.",
                         field_types=val_dict, images=img_paths, default_text='Advanced Settings', default_fields=advanced_dict,
                         checkboxes=check_list, misc_files=misc_paths, folders = folder_paths, scroll=True)

        self.current_step = App.SEG_IP

        self.add_page(seg_ip)

        self.del_page(1)

        self.pathway = App.FROM_SEG

        self.cont_button.config(command=self.next_step, state='normal')


    def start_metrics(self):
        print('time 4 metrics')