import tkinter as tk
from tkinter import filedialog as fd
import copy
from PIL import Image, ImageTk
import numpy as np
import cv2
from tkinter.font import Font, nametofont

"""
@author Tanya Roysam
roysamt@gatech.edu
Page.py
Child class of tk.Frame
Creates page for user input and output display.
* Float value inputs
* Image path inputs (currently only .png, .tif)
* Save paths 
"""


class Page(tk.Frame):
    page_total = 0

    def __init__(self, root, title, field_types=None, file_saves=None, images=None, default_text = None, default_fields = None,
                 checkboxes = None, misc_files = None, folders = None, scroll=False):
        """
        Constructor for page class. Descendant of tk.Frame

        :param root: Master window this frame is in. Navigator object, usually
        :param title: Title of the page
        :param field_types: Dictionary {name of field : type of field}. use 'str' or 'float'
        :param file_saves: types of files required for saves
        :param images: image file inputs
        :param default_text: text with default input fields
        :param default_fields: dict {field name : default field value} RIGHT NOW ALL DEFAULT FIELDS MUST BE FLOAT
        :param checkboxes: list of prompts for checkbox data
        """
        super().__init__(root)
        # increment 1 to total number of pages

        self.MAX_WIDTH = int(2*root.winfo_screenwidth()/3) - 200


        Page.page_total = Page.page_total + 1
        self.pgc = Page.page_total

        self.r = root

        # display title at top
        self.title = title

        title_font = Font(size=16)
        tk.Label(self,text=title, font=title_font).pack(side=tk.TOP)

        self.default_font = nametofont('TkDefaultFont')




        if scroll:
            self.scrollbar = tk.Scrollbar(self, orient=tk.VERTICAL)
            self.scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
            self.xscrollbar = tk.Scrollbar(self, orient=tk.HORIZONTAL)
            self.xscrollbar.pack(side=tk.BOTTOM, fill=tk.X)
            self.canvas = tk.Canvas(self, highlightthickness=0)
            self.scrollbar.config(command=self.canvas.yview)
            self.xscrollbar.config(command=self.canvas.xview)
            self.canvas.config(yscrollcommand=self.scrollbar.set, xscrollcommand=self.xscrollbar.set)
            self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

            self.canvas.bind('<Enter>', self.bind_mouse)
            self.canvas.bind('<Leave>', self.unbind_mouse)

            self.c_frame = tk.Frame(self.canvas)
            anchor_x = int(root.winfo_screenwidth() / 3)
            anchor_y = int(root.winfo_screenheight() / 3)



            self.canvas.create_window(anchor_x, anchor_y, window=self.c_frame, anchor=tk.CENTER)
            self.c_frame.bind('<Configure>', self.on_frame_configure)
        else:
            self.canvas = tk.Canvas(self, highlightthickness=0)
            self.c_frame = tk.Frame(self.canvas)
            self.canvas.pack()
            self.c_frame.pack()

        if field_types:
            # frame for entry fields
            self.field_frame = tk.LabelFrame(self.c_frame, text='Value Entry', padx=10, pady=10)
            self.field_frame.grid_columnconfigure(0, weight=1)
            self.field_frame.grid_columnconfigure(1, weight=1)

            self.field_frame.pack_propagate(0)
            # generating field dictionaries for storage & generating fields/labels
            self.field_types = field_types
            field_names = list(field_types.keys())
            self.field_dict = {field_names[i] : None for i in range(len(field_names))}
            self.field_value = {field_names[i] : None for i in range(len(field_names))}

            self.diag_text = tk.StringVar(self.field_frame, value=' ')
            self.fieldgen()


            self.field_frame.pack(pady=20, fill='both', expand='yes')

        if checkboxes:
            self.check_frame = tk.LabelFrame(self.c_frame, text='Options', padx=10, pady=10)

            #self.checkframe.pack_propagate(0)
            self.check_vars = {checkboxes[i]:None for i in range(len(checkboxes))}
            self.check_vals = {checkboxes[i]:False for i in range(len(checkboxes))}
            self.create_check()

        if file_saves:
            # file input frame
            self.file_frame = tk.LabelFrame(self.c_frame, text='Save Outputs', padx=10, pady=10)

            self.file_frame.pack_propagate(0)
            # file input dictionaries and generating file input dialogs/buttons
            self.file_names = {file_saves[i]:None for i in range(len(file_saves))}
            self.file_disps = {file_saves[i]:tk.StringVar(self.file_frame, value ='Selected file:') for i in range(len(file_saves))}
            self.longest_file = max(file_saves, key=len)
            self.create_savefile()
            self.file_frame.pack(pady=20, fill='both', expand='yes')

        if images:

            # image file input frame
            self.imgfile_frame = tk.LabelFrame(self.c_frame, text='Upload Images', padx=10, pady=10)

            self.imgfile_frame.pack_propagate(0)
            # image file input dictionaries and generating file input dialogs/buttons
            self.imgfile_names = {images[i]: None for i in range(len(images))}
            self.imgfile_disps = {images[i]: tk.StringVar(self.imgfile_frame, value='Selected file:') for i in range(len(images))}
            self.fileimg = {images[i]: None for i in range(len(images))}
            self.images = {images[i]: None for i in range(len(images))}
            self.longest_img = max(images, key=len)
            self.create_imagefile()
            self.imgfile_frame.pack(pady=20, fill='both', expand='yes')



        if misc_files:
            self.miscfile_frame = tk.LabelFrame(self.c_frame, text='Upload Other Files', padx=10, pady=10)

            self.miscfile_frame.pack_propagate(0)
            self.miscfile_names = {misc_files[i]: None for i in range(len(misc_files))}
            self.miscfile_disps = {misc_files[i]: tk.StringVar(self.miscfile_frame, value='Selected file:') for i in range(len(misc_files))}
            self.longest_misc = max(misc_files, key=len)
            self.create_miscfile()
            self.miscfile_frame.pack(pady=20, fill='both', expand='yes')

        if folders:
            self.folder_frame = tk.LabelFrame(self.c_frame, text='Select Output Folder', padx=10, pady=10)
            self.folder_frame.pack_propagate(0)

            self.folder_names = {folders[i]: None for i in range(len(folders))}
            self.folder_disps = {folders[i]: tk.StringVar(self.folder_frame, value='Selected file:') for i in range(len(folders))}
            self.longest_folder = max(folders, key=len)

            self.create_folder()
            self.folder_frame.pack(pady=20, fill='both', expand='yes')

        if default_fields:
            # Idea: Open Advanced Settings button
            self.default_vals = default_fields
            self.default_fields = {list(default_fields.keys())[i] : None for i in range(len(default_fields))}
            self.default_b = tk.Button(self.c_frame, text=default_text, command=self.default_open)
            self.default_b.pack()
            self.default_frame = tk.LabelFrame(self.c_frame, text='Settings', padx=10, pady=10)

            self.default_frame.pack_propagate(0)

            self.default_frame.grid_columnconfigure(0, weight=1)
            self.default_frame.grid_columnconfigure(1, weight=1)
            self.create_default()

            self.def_open = False


    # SCROLLBAR METHODS

    def on_frame_configure(self, event):
        self.canvas.config(scrollregion=self.canvas.bbox('all'))

        self.canvas.config(yscrollcommand=self.scrollbar.set, xscrollcommand=self.xscrollbar.set)
        self.scrollbar.config(command=self.canvas.yview)
        self.xscrollbar.config(command=self.canvas.xview)
    def bind_mouse(self, event):
        self.canvas.bind_all('<4>', self.on_mousewheel)
        self.canvas.bind_all('<5>', self.on_mousewheel)
        self.canvas.bind_all('<MouseWheel>', self.on_mousewheel)
    def unbind_mouse(self, event):
        self.canvas.unbind_all('<4>')
        self.canvas.unbind_all('<5>')
        self.canvas.unbind_all('<MouseWheel>')
    def on_mousewheel(self, event):
        if event.num == 4 or event.delta > 0:
            self.canvas.yview_scroll(-1, 'units')
        elif event.num == 5 or event.delta < 0:
            self.canvas.yview_scroll(1, 'units')

    def shorten_path(self, path, mlen):
        sflen = self.default_font.measure('Selected file: ...')
        mlen -= sflen
        while self.default_font.measure(path) > mlen:
            ind = path.find('/')
            if ind == -1: #return os.path.basename(path)
                while self.default_font.measure(path) > mlen:
                    path = path[1:]
                return path
            if len(path) > ind+1:
                path = path[ind+1:]
            else:
                return path[-20:]
        return path




    def default_open(self):
        if not self.def_open:
            self.default_frame.pack(pady=20, fill='both', expand='yes')
            self.def_open = True
        else:
            self.default_frame.pack_forget()
            self.def_open = False


    def nav(self):
        """
        'Navigate' to this page

        :return: (int) number of page
        """
        self.lift()

    def fieldgen(self):
        """
        Generate input fields. ONLY call in constructor

        :return: None
        """
        row = 1

        # iterate through fields to generate
        for k in self.field_dict:

            # field label to left of entry
            tk.Label(self.field_frame, text=k).grid(row=row, column=0, sticky=tk.E)

            # create and position entry
            en = tk.Entry(self.field_frame)
            en.grid(row=row, column=1, sticky=tk.W)

            # update field_dict to reflect new entry field, increment row
            self.field_dict.update({k : en})
            row = row + 1

        # Create save button- will collect all value inputs and check for validity
        tk.Button(self.field_frame, text ='Save', command = self.get_field).grid(row=row, columnspan=2, sticky=tk.N)
        tk.Label(self.field_frame, textvariable=self.diag_text).grid(row=row + 1, columnspan=2, sticky=tk.N)

    def get_field(self):
        """
        Collect all inputs from Entry fields & cast, save

        :return: None
        """
        # iterate through all Entry objects saved in dictionary
        for k, val in self.field_dict.items():
            v = val.get()

            # if there is input in the field
            if v != '':
                if self.field_types[k] == 'float':
                    # try to cast as float if type is specified float. if not, display incorrect datatype
                    try:
                        v = float(v)
                        self.field_value.update({k: v})
                        self.diag_text.set(' ')
                    except ValueError:
                        self.diag_text.set('incorrect datatype!')

                else:
                    # if it's not a float, it's a string, so it can be updated without casting
                    self.field_value.update({k: v})
            else:
                # input being '' indicates that no input was made
                self.diag_text.set('please enter value to continue.')
            # update the labels
            self.update_idletasks()

    def create_savefile(self):
        """
        Generate file save fields. ONLY run in constructor

        :return: None
        """
        row = 0
        # iterate through all the files needed
        for label in self.file_names.keys():
            # Need to define this lambda function individually for every button
            # because variable label helps identify which button has been pressed
            def button_cmd(lbl=label):
                self.get_savefile(lbl)
            # create button with command opening up file dialogue and label confirming file input
            tk.Button(self.file_frame, text='Save ' + label, command=button_cmd).grid(row=row, column=0, sticky=tk.W)
            tk.Label(self.file_frame, textvariable=self.file_disps[label]).grid(row=row, column=1, sticky=tk.W)
            row = row + 1

    def create_folder(self):
        """
        Generate folder save fields. ONLY run in constructor

        :return: None
        """
        row = 0
        # iterate through all the files needed
        for label in self.folder_names.keys():
            # Need to define this lambda function individually for every button
            # because variable label helps identify which button has been pressed
            def button_cmd(lbl=label):
                self.dir_save(lbl)
            # create button with command opening up file dialogue and label confirming file input
            tk.Button(self.folder_frame, text='Select ' + label, command=button_cmd).grid(row=row, column=0, sticky=tk.W)
            tk.Label(self.folder_frame, textvariable=self.folder_disps[label]).grid(row=row, column=1, sticky=tk.W)
            row = row + 1


    def get_savefile(self, file):
        """
        Runs when button is clicked. Retrieve file path and save

        :param file: (string) Identifier of which file's path is being retrieved. File path retrieved is saved into
        dictionary attribute under the key corresponding to this parameter
        :return: None
        """
        # open file dialog
        fin = fd.asksaveasfilename()
        # update filenames dict to = fin
        self.file_names.update({file:fin})

        # change respective label
        max_len = self.MAX_WIDTH - self.default_font.measure(self.longest_file)

        if fin != None:
            if self.default_font.measure(fin) > max_len:
                flen = self.shorten_path(fin, max_len)
                self.file_disps[file].set('Selected file: ...' + flen)
            else:
                self.file_disps[file].set('Selected file: ' + fin)
        # update labels
        self.update_idletasks()

    def dir_save(self, file):
        """
        open filedialog to save directory

        :param file: (string) Identifier of which file's path is being retrieved. File path retrieved is saved into
        dictionary attribute under the key corresponding to this parameter
        :return: None
        """
        fin = fd.askdirectory(mustexist=True)
        fin += '/'
        # update filenames dict to = fin
        self.folder_names.update({file: fin})
        max_len = self.MAX_WIDTH - self.default_font.measure(self.longest_folder)

        # change respective label
        if fin != '/':
            if self.default_font.measure(fin) > max_len:
                flen = self.shorten_path(fin, max_len)
                self.folder_disps[file].set('Selected file: ...' + flen)
            else:
                self.folder_disps[file].set('Selected file: ' + fin)
        else:
            self.folder_names.update({file:None})
        # update labels
        self.update_idletasks()

    def get_info(self):
        """
        Getter for user input. Returns deep copies
        :return: (dicts) field values, filenames, imgfilenames, default values, checkbox values
        """
        try:
            fval = copy.deepcopy(self.field_value)
        except AttributeError:
            fval = {}
        try:
            fn = copy.deepcopy(self.file_names)
        except AttributeError:
            fn = {}
        try:
            imgfn = copy.deepcopy(self.imgfile_names)
        except AttributeError:
            imgfn = {}
        try:
            defval = copy.deepcopy(self.default_vals)
        except AttributeError:
            defval = {}
        try:
            checkval = copy.deepcopy(self.check_vals)
        except AttributeError:
            checkval = {}
        try:
            miscfn = copy.deepcopy(self.miscfile_names)
        except AttributeError:
            miscfn = {}
        try:
            foldfn = copy.deepcopy(self.folder_names)
        except AttributeError:
            foldfn = {}

        return fval, fn, imgfn, defval, checkval, miscfn, foldfn

    def create_imagefile(self):
        """
        Generate image input fields and accompanying thumbnails of selected image

        :return: None
        """
        row = 0

        # iterate through all the files needed
        for label in self.imgfile_names.keys():

            # Need to define this lambda function individually for every button
            # because variable label helps identify which button has been pressed
            def button_cmd(lbl=label):
                """
                this is a lambda. local
                :param lbl: which button pressed
                :return: none
                """
                self.get_imagefile(lbl)

            # create button with command opening up file dialogue and label confirming file input
            tk.Button(self.imgfile_frame, text='Upload ' + label, command=button_cmd).grid(row=row * 2, column=0, sticky=tk.W)

            # Label for what file you opened
            tk.Label(self.imgfile_frame, textvariable=self.imgfile_disps[label]).grid(row=row * 2, column=1, sticky=tk.W)
            k = tk.Frame(self.imgfile_frame)
            k.grid(row=row * 2 + 1,column=1)
            self.fileimg.update({label: k})
            row = row + 1

    def get_imagefile(self, file):
        """
        Open filedialog input, load image for thumbnail

        :param file:
        :return:
        """
        # open file dialog
        fin = fd.askopenfilename()

        # update filenames dict to = fin
        self.imgfile_names.update({file: fin})

        max_len = self.MAX_WIDTH - self.default_font.measure(self.longest_img)
        # change respective label
        if fin != None:
            if self.default_font.measure(fin) > max_len:
                flen = self.shorten_path(fin, max_len)
                self.imgfile_disps[file].set('Selected file: ...' + flen)
            else:
                self.imgfile_disps[file].set('Selected file: ' + fin)

        # check if it's a compatible file type
        if (fin[-4:] != ".tif") & (fin[-4:] != ".png"):
            # it is not
            self.imgfile_names.update({file:None})
            self.imgfile_disps[file].set('Selected file: That is not a compatible filetype. Please try again.')
        else:
            # it is valid. destroy any existing image thumbnail and get raw image
            for w in self.fileimg[file].winfo_children():
                w.destroy()
            rawimg = cv2.imread(fin,0)

            # convert raw image to PIL Image, resize to thumbnail max size, and convert to PIL ImageTk
            rawimg= Image.fromarray(
                    (cv2.normalize(rawimg, None, alpha=0, beta=255, norm_type=cv2.NORM_MINMAX)).astype(np.uint8))
            rw, rh = rawimg.size
            if rw > rh:
                scalar = rw/200
                nw = 200
                nh = int(rh / scalar)
            else:
                scalar = rh/200
                nh = 200
                nw = int(rw / scalar)
            rawimg = rawimg.convert(mode="L")
            rawimg = rawimg.resize((nw, nh))
            img = ImageTk.PhotoImage(rawimg)

            # update dictionary attribute so image renders, and display
            self.images.update({file:img})
            tk.Label(self.fileimg[file], image=self.images[file]).pack()

        self.update_idletasks()

    def check_input(self):
        """
        Check if all input is present
        If input is in dictionaries, it is valid input bc above functions check for validity

        :return: (bool) if all inputs were put in
        """
        complete = True
        try:
            for i in self.field_value.values():
                if i==None: return False
        except AttributeError:
            complete=True
        try:
            for i in self.file_names.values():
                if i == None: return False
        except AttributeError:
            complete = True
        try:
            for i in self.imgfile_names.values():
                if i == None: return False
        except AttributeError:
            complete = True
        try:
            for i in self.default_vals.values():
                if i == None: return False
        except AttributeError:
            complete = True
        try:
            for i in self.miscfile_names.values():
                if i == None: return False
        except AttributeError:
            complete = True
        try:
            for i in self.folder_names.values():
                if i == None: return False
        except AttributeError:
            complete = True
        return complete

    def create_default(self):
        """
        generate input fields with default settings

        :return: None
        """
        # iterate through fields to generate
        for ind, (k,val) in enumerate(self.default_vals.items()):
            # field label to left of entry
            tk.Label(self.default_frame, text=k).grid(row=ind + 1, column=0, sticky=tk.E)

            # create and position entry
            en = tk.Entry(self.default_frame)
            en.insert(tk.END, str(val))
            en.grid(row=ind+1, column=1, sticky=tk.W)

            # update field_dict to reflect new entry field, increment row
            self.default_fields.update({k: en})

        # Create save button- will collect all value inputs and check for validity
        tk.Button(self.default_frame, text='Save', command=self.get_default).grid(row=ind + 2, columnspan=2, sticky=tk.N)

    def get_default(self):
        """
        fetch inputs in default fields

        :return: None
        """
        # iterate through all Entry objects saved in dictionary
        for k, val in self.default_fields.items():
            v = val.get()


            # if there is input in the field
            if v != '':
                try:
                    v = float(v)
                    if v.is_integer(): v=int(v)
                    self.default_vals.update({k: v})
                    self.diag_text.set(' ')
                except ValueError:
                    self.diag_text.set('incorrect data type!')
            else:
                # input being '' indicates that no input was made
                self.diag_text.set('please enter value to continue.')
            # update the labels
            self.update_idletasks()

    def create_check(self):
        """
        generate checkbox fields

        :return: None
        """
        for k, val in self.check_vars.items():
            var = tk.IntVar(value=0)
            tk.Checkbutton(self.check_frame, text=k, variable=var, command=self.get_check, onvalue=1, offvalue=0).pack()
            self.check_vars.update({k:var})
        self.check_frame.pack(pady=20, fill='both', expand='yes')
        self.update_idletasks()


    def get_check(self):
        """
        retrieve checkbox data

        :return: None
        """

        for k, var in self.check_vars.items():
            self.check_vals.update({k: bool(var.get())})

    def create_miscfile(self):
        """
        Generate generic file input fields

        :return: None
        """
        row = 0

        # iterate through all the files needed
        for label in self.miscfile_names.keys():

            # Need to define this lambda function individually for every button
            # because variable label helps identify which button has been pressed
            def button_cmd(lbl=label):
                """
                this is a lambda. local
                :param lbl: which button pressed
                :return: none
                """
                self.get_miscfile(lbl)

            # create button with command opening up file dialogue and label confirming file input
            tk.Button(self.miscfile_frame, text='Upload ' + label, command=button_cmd).grid(row=row * 2, column=0, sticky=tk.W)

            # Label for what file you opened
            tk.Label(self.miscfile_frame, textvariable=self.miscfile_disps[label]).grid(row=row * 2, column=1, sticky=tk.W)
            row = row + 1

    def get_miscfile(self, file):
        """
        Open filedialog input

        :param file:
        :return:
        """
        # open file dialog
        fin = fd.askopenfilename()

        # update filenames dict to = fin
        self.miscfile_names.update({file: fin})
        max_len = self.MAX_WIDTH - self.default_font.measure(self.longest_misc)

        # change respective label
        if fin != None:
            if self.default_font.measure(fin) > max_len:
                flen = self.shorten_path(fin, max_len)
                self.miscfile_disps[file].set('Selected file: ...' + flen)
            else:
                self.miscfile_disps[file].set('Selected file: ' + fin)

        # check if it's a compatible file type
        if (fin[len(fin)-4:].lower() != ".csv") & (fin[len(fin)-4:].lower() != ".txt") & (fin[len(fin)-6:].lower() != ".imzml"): # TODO: add any other file types
            # it is not
            self.miscfile_names.update({file:None})
            self.miscfile_disps[file].set('Selected file: That is not a compatible filetype. Please try again.')
            #Comma-separated list of ions (e.g. 883.576,885.573,887.57)


        self.update_idletasks()