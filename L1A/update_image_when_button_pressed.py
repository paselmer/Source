from tkinter import *
from tkinter import ttk
from PIL import ImageTk, Image

# NOTES: You need to make a python list of "PIL Photo objects."
#        When you click the button and enter the function, you
#        open the updated image as a "photo object," overwrite
#        the element of the list, then update the label.

def change_img(*args):
    print('hello')
    photos[1] = ImageTk.PhotoImage(Image.open(im2))
    label.config(image=photos[1])
    
root = Tk()
root.title("Experiments")

mainframe = ttk.Frame(root, width=1350, height=350, borderwidth= 2, relief='sunken')
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
mainframe.columnconfigure(0, weight=1)
mainframe.rowconfigure(0, weight=1)

im1 = '/home/cpl/camal/source/def.gif'
im2 = '/home/cpl/camal/source/new_image.png'
photos = [ ImageTk.PhotoImage(Image.open(im1)), None ]
label = ttk.Label(mainframe, image=photos[0])
label.grid(column=1, row=0, sticky=E)

b = Button(mainframe,text='load the image',command=change_img).grid(column=0, row=0, sticky=N)

# The following commands tell the "mainframe" how to change the shape of its
# columns if the user resizes the window.
for child in mainframe.winfo_children(): child.grid_configure(padx=5,pady=5)


root.mainloop()
