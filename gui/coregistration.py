import tkinter as tk
import App


def main():
    master = tk.Tk()
    master.title('Coregistration')

    master.geometry('%sx%s' % (int(2*master.winfo_screenwidth()/3), master.winfo_screenheight()-150))

    app = App.App(master)
    app.mainloop()


if __name__ == "__main__":
    main()
