import tkinter as tk
import Output_Page as pf
import threading
import queue
import logging

"""
@author Tanya Roysam
roysamt@gatech.edu
Navigator.py
Child class of tk.Frame. Overall navigation object that contains Page objects.
Handles function of overall GUI. 
* Navigation from page to page
* Checking and fetching user input
* Initiating secondary threads and handling output of threads
"""

class Navigator(tk.Frame):
    """
    Child class of tk.Frame. Overall navigation object that contains Page objects.
    Handles function of overall GUI.
    * Navigation from page to page
    * Checking and fetching user input
    * Initiating secondary threads and handling output of threads
    """

    def __init__(self, root, pagelist, methods, output_titles, argnames):
        """
        Constructor for overall app

        :param root: (tk) Master window
        :param pagelist: (list, Page) list of initial pages to add
        """

        # parent constructor
        super(Navigator, self).__init__(root)

        self.r = root
        self.pack()

        # reference to root object, page marker, list of pages
        self.current_page = 1
        self.page_list = pagelist

        self.MAX_HEIGHT = root.winfo_screenheight() - 150
        self.MAX_WIDTH = int(2*root.winfo_screenwidth() / 3)

        page_height = int(self.MAX_HEIGHT * 7 / 8)
        button_height = int(self.MAX_HEIGHT / 8)

        # container frames for pages, navigation buttons
        self.page_container = tk.Frame(self, width = self.MAX_WIDTH, height = page_height)
        self.button_frame = tk.Frame(self, width = self.MAX_WIDTH, height = button_height)


        # place pages into container frame
        for pg in self.page_list:
            pg.place(in_=self.page_container, relheight=1, relwidth=1)

        # navigate to first page
        self.page_list[0].nav()
        self.page_container.pack()

        # Add in page navigation buttons
        self.next_b = tk.Button(self.button_frame, text='Next Page', command=self.next_page)
        self.next_b.grid(row=1, column=2)
        self.prev_b = tk.Button(self.button_frame, text='Previous Page', command=self.back_page)
        self.prev_b.grid(row=1, column=1)
        self.button_frame.pack()

        if len(self.page_list) == 1:
            self.next_b.config(state='disabled')
            self.prev_b.config(state='disabled')

        # Threading object
        self.thread = None

        # Output queue (shared state variable)
        self.out = queue.Queue()

        # Output attribute. fetched from self.out
        self.current_out = None

        # Shared state integer variable
        self.prog_var = tk.IntVar()

        # Current loading page
        self.lp = None

        # continue button at bottom of page
        self.cont_button = tk.Button(self.button_frame, text='Continue', state='disabled')
        self.cont_button.grid(row=2, column=1, columnspan=2, sticky=tk.N)
        self.cont_text = tk.StringVar(self.button_frame, value=' ')
        tk.Label(self, textvariable=self.cont_text).pack(side=tk.BOTTOM)

        self.pack()

        self.methods = methods
        self.current_step = 1

        # USER CHANGE
        self.output_titles = output_titles

        self.inputs = {'pv':self.prog_var}

        self.arg_names = argnames

        self.failed = tk.IntVar()
        self.exs = queue.Queue()


    def check(self, ip_nums):
        """
        Check inputs from input pages for validity

        :param ip_nums: Page(s) to check for inputs
        :return: (bool) True if all inputs are valid, False if not
        """
        # ip nums starts from 1, input num starts from 1
        yes = True

        for i in ip_nums:
            p = self.page_list[i - 1]
            yes = yes & p.check_input()

        if yes:
            for i in ip_nums:
                p = self.page_list[i - 1]
                # output dicts vals, saves, imgpaths, defvals, checkvals, miscpaths
                out_dicts = p.get_info()
                for dict in out_dicts:
                    self.inputs.update(dict)

        return yes

    def get_args(self, input_num):
        """
        Get args from a given page

        :param input_num: (int) page number to retrieve data from
        :return: (list) args
        """
        args = []

        argname_list = self.arg_names[input_num - 1]

        for name in argname_list:
            args.append(self.inputs[name])

        return args


    def add_page(self, pg_new):
        """
        add page (Page or Output_Page)

        :param pg_new: (Page) new page to be added
        :return: None
        """
        # Append to pagelist, place within container frame

        self.page_list.append(pg_new)
        pg_new.place(in_=self.page_container, relheight=1, relwidth=1)

        # Set current page and navigate to page
        self.current_page = len(self.page_list)

        if len(self.page_list) == 1:
            self.prev_b.config(state='disabled')
            self.next_b.config(state='disabled')
        else:
            self.prev_b.config(state='normal')



        pg_new.nav()


        # update
        self.update_idletasks()

    def del_page(self, num):
        """
        Delete page (specified by number, starting at 1)

        :param num: (int) page number to delete
        :return: None
        """

        # return page to be deleted and destroy it, remove from pagelist
        dp = self.page_list[num - 1]
        dp.destroy()
        del self.page_list[num - 1]

        # navigate to page right before it
        if self.current_page==num:
            if num > 1:
                self.current_page = num - 1
                self.page_list[num - 2].nav()

        if num == 1:
            self.current_page = num
            self.page_list[num - 1].nav()

        if self.current_page == 1:
            self.prev_b.config(state='disabled')
            if len(self.page_list) > 1:
                self.next_b.config(state='normal')

        # update
        self.update_idletasks()


    def next_page(self):
        """
        Navigate to next page

        :return: None
        """

        if self.current_page < len(self.page_list):
            self.current_page = self.current_page + 1
            self.page_list[self.current_page - 1].nav()

        if self.current_page == len(self.page_list):
            self.next_b.config(state='disabled')

        if self.current_page > 1:
            self.prev_b.config(state='normal')

    def back_page(self):
        """
        Navigate to previous page

        :return: None
        """

        if self.current_page > 1:
            self.current_page = self.current_page - 1
            self.page_list[self.current_page - 1].nav()

        if self.current_page == 1:
            self.prev_b.config(state='disabled')

        if self.current_page < len(self.page_list):
            self.next_b.config(state='normal')



    def loading_page(self, name, iters, det):
        """
        Create and add loading page

        :param name: (string) name of loading page (usually task being executed)
        :param iters: (int) number of steps needed (for progress bar)
        :return: None
        """
        self.prog_var.set(0)
        self.lp = pf.Output_Page(self, name, '0', [iters, self.prog_var, det])
        self.add_page(self.lp)

    def del_loading(self):
        """
        Remove current loading page

        :return: None
        """
        if self.lp is not None:
            self.del_page(len(self.page_list))
            self.lp = None

    def output(self, num, method = None):
        # Step numbers start at 1
        self.add_page(pf.Output_Page(self, self.output_titles[num - 1], num, self.current_out, method))


    def check_status(self, nm):
        """
        Check status of thread.
        If it is dead (aka finished executing) save output from Queue, reconfigure button, generate output page.

        :param nm: (list) which step of pipeline, method for checkbox option (None if no checkbox)
        :return: None
        """
        # If the thread is alive, keep calling checkstatus every 250 ms and update loading page (for prog bar)
        if self.thread.is_alive():
            self.after(400, self.check_status, nm)
            if self.lp:
                self.lp.update_idletasks()

        else:
            # Depending on which step, get output from queue and generate page with output
            # Then configure continue button to point to next thread
            if self.failed.get()==1:
                self.error_handling(self.exs.get())
            else:
                num = nm[0]
                method = nm[1]
                out = self.out.get()
                self.del_loading()
                self.current_out = out

                self.output(num, method)

                # Config contbutton somewhere
                self.cont_button.config(state='normal')

            # update
            self.update_idletasks()

    def error_handling(self, ex):
        """
        Handle error by shutting down program

        :param ex: (Exception) exception thrown
        :return: None
        """
        for w in self.page_container.winfo_children():
            w.destroy()
        self.page_list = []
        self.cont_button.config(state='disabled')
        self.next_b.config(state='disabled')
        self.prev_b.config(state='disabled')
        pg = pf.Output_Page(self.r, 'Error', -1, [ex])
        self.add_page(pg)

    def create_thread(self, args, method, pvnum, lptitle, mnum, method_cbox=None):
        """
        Create and start second thread for method to run in

        :param args: (list) args for function step
        :param method: (method) method to be run in thread
        :param pvnum: (int) number of total steps for progress variable
        :param lptitle: (string) title of loading page
        :param mnum: (int) method number
        :param method_cbox: (method) method associated with checkbox
        :return:
        """
        def run_thread():
            # Local. acts like lambda
            # needed for error handling
            try:
                self.out.put(method(*args))
            except Exception as ex:
                self.failed.set(1)
                self.exs.put(ex)
                logging.exception('Program was shut down due to the following exception:')
            finally:
                pass

        self.thread = threading.Thread(target=run_thread)

        if pvnum == 0:
            d = False
        else:
            d = True

        if lptitle is not None:
            self.loading_page(lptitle, pvnum, det=d)

        self.cont_button.config(state='disabled')

        self.thread.start()

        self.check_status([mnum, method_cbox])
