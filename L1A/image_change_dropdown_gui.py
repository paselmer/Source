from tkinter import *
from tkinter import ttk

def change_img(*args):
    print('hello')
    label.config(image=photos[int(imFiles.get())])
    
root = Tk()
root.title("Experiments")

mainframe = ttk.Frame(root, width=1350, height=350, borderwidth= 2, relief='sunken')
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
mainframe.columnconfigure(0, weight=1)
mainframe.rowconfigure(0, weight=1)

imFiles = StringVar()
imFiles.set('Image')
files=['0','1'] # Number corresponds to list index
imDropDown = OptionMenu(mainframe, imFiles, *files)
imDropDown.config(font=("Times",16,"italic"))
imDropDown["menu"].config(font=("Times",16,"italic"))
imDropDown.grid(column=0,row=0)

im1 = '/home/cpl/camal/source/new_image.gif'
im2 = '/home/cpl/camal/source/def.gif'
photos = (PhotoImage(file=im1),PhotoImage(file=im2))
label = ttk.Label(mainframe, image=photos[0])
label.grid(column=1, row=1, sticky=E)

imFiles.trace("w",change_img)

# The following commands tell the "mainframe" how to change the shape of its
# columns if the user resizes the window.
for child in mainframe.winfo_children(): child.grid_configure(padx=5,pady=5)


root.mainloop()
