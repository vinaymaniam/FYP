import Tkinter as tk
import ttk, tkFileDialog, tkMessageBox
from PIL import Image, ImageTk
import subprocess
import math
import time
import os
import traceback
from datetime import timedelta
from sys import platform, exc_info

class ScrollableImage(tk.Frame):
    def __init__(self, parent, **kwargs):
        # normalization factor for mouse scrolling
        if platform=="darwin":
            self.scroll_norm = 1 
        else:
            self.scroll_norm = 120
        # initialization of Label embedding a scrollable canvas
        tk.Frame.__init__(self,parent,**kwargs)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.grid(row=0, column=0, sticky=tk.W+tk.S+tk.E+tk.N)
        self.canvas = tk.Canvas(self, highlightthickness=0)
        self.canvas.grid()
        self.image_canvas = self.canvas.create_image(0, 0, anchor=tk.N+tk.W)
        self.xbar = tk.Scrollbar(self, orient=tk.HORIZONTAL)
        self.xbar.grid(row=1, column=0, sticky=tk.E+tk.W)
        self.ybar = tk.Scrollbar(self, orient=tk.VERTICAL)
        self.ybar.grid(row=0, column=1, sticky=tk.N+tk.S)
        self.canvas.bind("<Enter>", self.on_entry)
        self.canvas.bind("<Leave>", self.on_exit)
    
    def set_image(self, data):
        self.tk_data = ImageTk.PhotoImage(data)
        self.canvas.itemconfig(self.image_canvas, image=self.tk_data)
        self.canvas.config(height=self.tk_data.height(), width=self.tk_data.width())
        self.canvas.config(height=self.tk_data.height(), width=self.tk_data.width(), 
        xscrollcommand=self.xbar.set, yscrollcommand=self.ybar.set, 
            scrollregion=self.canvas.bbox(tk.ALL))
        self.xbar.config(command=self.canvas.xview)
        self.ybar.config(command=self.canvas.yview)
    
    def get_canvas(self):
        return self.canvas
        
    def on_vertical_scroll(self, event):
        self.canvas.yview_scroll(int(-1*(event.delta/self.scroll_norm)), "units")
    
    def on_horizontal_scroll(self, event):
        self.canvas.xview_scroll(int(-1*(event.delta/self.scroll_norm)), "units")
    
    def on_entry(self, event):
        if platform=="linux" or platform=="linux2":
            self.canvas.bind_all("<Button-4>", self.on_vertical_scroll)
            self.canvas.bind_all("<Button-5>", self.on_horizontal_scroll)
        else:
            self.canvas.bind_all("<MouseWheel>", self.on_vertical_scroll)
            self.canvas.bind_all("<Shift-MouseWheel>", self.on_horizontal_scroll)
        
    def on_exit(self, widget):
        if platform=="linux" or platform=="linux2":
            self.canvas.unbind_all("<Button-4>")
            self.canvas.unbind_all("<Button-5>")
        else:
            self.canvas.unbind_all("<MouseWheel>")
            self.canvas.unbind_all("<Shift-MouseWheel>")
        

class FreshGui(tk.Frame):
    def __init__(self, root, *args, **kwargs):
        # internal variables
        self.x = 0
        self.y = 0
        self.block_size = 64
        self.level = 1
        self.diag_reg = 1
        self.lin_map = 1
        self.test_mode = 0
        self.fast_profile = 0
        
        self.title = "FRESH Super-Resolution"
        self.warning_8X_not_shown = True
        
        # normalization factor for mouse scrolling
        if platform=="darwin":
            self.scroll_norm = 1 
        else:
            self.scroll_norm = 120 # 120 for win and linux
        
        # initialize root container
        tk.Frame.__init__(self, root, *args, **kwargs)
        self.root = root
        self.root.title("FRESH Single-Image Super-Resolution")
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        
        # paned window for top row
        self.paned_window = tk.PanedWindow(self.root, sashpad=5, sashwidth=4, 
            sashrelief=tk.GROOVE, showhandle=True, opaqueresize=False)
        self.paned_window.grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky=tk.W+tk.S+tk.E+tk.N)
        
        # image container
        self.scroll_img_input = ScrollableImage(self.root)
        self.paned_window.add(self.scroll_img_input, stretch="always")
        # bind events on image
        self.rect = self.scroll_img_input.get_canvas().create_rectangle(0,0,0,0, width=0)
        self.prev_rect = self.scroll_img_input.get_canvas().create_rectangle(0,0,0,0, width=0)
        self.scroll_img_input.get_canvas().bind("<Enter>", self.on_entry_input_image)
        self.scroll_img_input.get_canvas().bind("<Leave>", self.on_exit_input_image)
        self.scroll_img_input.get_canvas().bind("<Motion>", self.on_move_input_image)
        self.scroll_img_input.get_canvas().bind("<Button-1>", self.on_click_input_image)
        
        # output block container
        self.notebook = ttk.Notebook(self.root)
        self.notebook.columnconfigure(0, weight=1)
        self.notebook.rowconfigure(0, weight=1)
        self.notebook.grid(row=0, column=0, sticky=tk.N+tk.S+tk.W+tk.E)
        self.paned_window.add(self.notebook, stretch="always")
        self.scroll_img_near = ScrollableImage(self.notebook)
        self.scroll_img_bic = ScrollableImage(self.notebook)
        self.scroll_img_fresh = ScrollableImage(self.notebook)
        self.notebook.add(self.scroll_img_fresh, text="FRESH Super-Resolution")
        self.notebook.add(self.scroll_img_bic, text="Bicubic Interpolation")
        self.notebook.add(self.scroll_img_near, text="Nearest Neighbour")
        
        # parameters container
        self.label_frame = tk.LabelFrame(self.root, text="Parameters")
        self.label_frame.grid(row=1, column=0, 
            padx=5, pady=5, ipadx=5, ipady=5, sticky=tk.W+tk.E)
        # select boxes
        tk.Label(self.label_frame, text="Block Size").grid(row=0, column=0, padx=5, sticky=tk.W)
        self.combo_block_size = ttk.Combobox(self.label_frame, state="readonly", 
            values=("32x32", "64x64", "128x128", "256x256"), width="12")
        self.combo_block_size.grid(row=1, column=0, padx=5, sticky=tk.W)
        tk.Label(self.label_frame, text="Resize Factor").grid(row=0, column=1, padx=5, sticky=tk.W)
        self.combo_level = ttk.Combobox(self.label_frame, state="readonly", 
            values=("2X", "4X", "8X"), width="12")
        self.combo_level.grid(row=1, column=1, padx=5, sticky=tk.W)
        tk.Label(self.label_frame, text="Profile").grid(row=0, column=2, padx=5, sticky=tk.W)
        self.combo_profile = ttk.Combobox(self.label_frame, state="readonly", 
            values=("Fastest", "Best Quality"), width="12")
        self.combo_profile.grid(row=1, column=2, padx=5, sticky=tk.W)
        # default values
        self.combo_block_size.current(1)
        self.combo_level.current(0)
        self.combo_profile.current(1)
        # bind events on control
        self.combo_block_size.bind("<<ComboboxSelected>>", self.change_block_size)
        self.combo_level.bind("<<ComboboxSelected>>", self.change_level)
        self.combo_profile.bind("<<ComboboxSelected>>", self.change_profile)
        
        # button container
        self.frame = tk.Frame(self.root)
        self.frame.columnconfigure(0, weight=1)
        self.frame.rowconfigure(0, weight=1)
        self.frame.grid(row=1, column=1, padx=5, pady=5, sticky=tk.W+tk.S+tk.E)
        # buttons
        self.button_save = tk.Button(self.frame, state=tk.DISABLED, text="Save Result", command=self.save_image)
        self.button_save.grid(row=0, column=0, padx=5, sticky=tk.W+tk.E)
        self.button_load = tk.Button(self.frame, text="Load New Image", command=self.load_image)
        self.button_load.grid(row=1, column=0, padx=5, sticky=tk.W+tk.E)
        self.button_upsample = tk.Button(self.frame, text="Upsample Image", command=self.upsample_image)
        self.button_upsample.grid(row=2, column=0, padx=5, sticky=tk.W+tk.E)
        
    def change_block_size(self, event):
        self.block_size = int(self.combo_block_size.get()[:self.combo_block_size.get().find("x")])
        self.combo_block_size.selection_clear()
        
    def change_level(self, event):
        self.level = int(math.log(float(self.combo_level.get()[:self.combo_level.get().find("x")]),2))
        if self.level==3 and self.warning_8X_not_shown:
            tkMessageBox.showwarning(self.title, "Upsampling factor 8X is experimental and might take a long time.")
            self.warning_8X_not_shown = False
        self.combo_level.selection_clear()
        
    def change_profile(self, event):
        if self.combo_profile.get()=="Fastest":
            self.diag_reg = 1
            self.lin_map = 0
            self.fast_profile = 1
        else:
            self.diag_reg = 1
            self.lin_map = 1
            self.fast_profile = 0
        self.combo_profile.selection_clear()
    
    def on_entry_input_image(self, event):
        # show outline of block to provess on mouse hovering
        self.scroll_img_input.get_canvas().itemconfig(self.rect, width=2)
    
    def on_exit_input_image(self, event):
        # set width to 0 (i.e., hide) block when mouse leaves
        self.scroll_img_input.get_canvas().itemconfig(self.rect, width=0)
    
    def get_coordinates(self, event):
        # compensate for offset due to scrollbars
        self.x = self.scroll_img_input.get_canvas().canvasx(event.x)
        self.y = self.scroll_img_input.get_canvas().canvasy(event.y)
    
    def overlay_block(self, rect_id, colour="yellow", line_width=2):
        # if block goes out of bounds then let's colour it red, otherwise it's yellow
        if self.x<self.block_size or self.y<self.block_size:
            colour = "red"
        # update position of rectangle
        self.scroll_img_input.get_canvas().itemconfig(rect_id, outline=colour, width=line_width)
        self.scroll_img_input.get_canvas().coords(rect_id, self.x, self.y, self.x-self.block_size, self.y-self.block_size)
    
    def on_move_input_image(self, event):
        self.get_coordinates(event)
        self.overlay_block(self.rect)
        
    def resource_path(self, relative_path):
        try:
            # PyInstaller creates a temp folder and stores path in _MEIPASS
            base_path = sys._MEIPASS
        except Exception:
            base_path = os.path.abspath(".")
        return os.path.join(base_path, relative_path)
    
    def upsample_data(self, block, verbose=True):
        if platform=="win32":
            binary = self.resource_path("fresh") + ".exe"
        else:
            binary = self.resource_path("fresh")
        if not os.path.isfile(binary):
            tkMessageBox.showerror(self.title, "FRESH binary not found.")
            raise RuntimeError("FRESH binary '" + binary + "' not found")
        # use temporary BMP files as input and output of FRESH binary
        tmp_filename_in = self.resource_path("input.bmp")
        tmp_filename_out = self.resource_path("output.bmp")
        block.convert("RGB").save(tmp_filename_in)
        # run FRESH
        cmd = "{} {} {} {} {} {} {} {}".format(binary, tmp_filename_in, tmp_filename_out, self.level, 
            self.diag_reg, self.lin_map, self.test_mode, self.fast_profile)
        if verbose:
            print("Upsampling: " + str(block.size) + "->" + str(tuple(i*(1<<self.level) for i in block.size)))
            print("Coordinates: " + str(tuple([int(self.x-self.block_size+1), int(self.y-self.block_size+1)])))
            print(cmd)
        start_time = time.time()
        try:
            fresh_binary_output = subprocess.check_output(cmd.split(), stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            tkMessageBox.showerror(self.title, e.output)
            raise RuntimeError("An error occurred during execution of FRESH")
        elapsed = (time.time() - start_time)
        if verbose:
            print(fresh_binary_output)
            print("Elapsed time: " + str(timedelta(seconds=elapsed)))
        # fetch FRESH results (open Pillow temporary and then copy 
        # becuse we want to retain the object but delete the file)
        tmp = Image.open(tmp_filename_out)
        self.block_fresh = tmp.copy()
        # standard upsampling methods
        self.block_bic = block.resize(self.block_fresh.size, Image.BICUBIC)
        self.block_near = block.resize(self.block_fresh.size, Image.NEAREST)
        # show upsampled blocks
        self.scroll_img_near.set_image(self.block_near)
        self.scroll_img_bic.set_image(self.block_bic)
        self.scroll_img_fresh.set_image(self.block_fresh)
        #enable save button
        self.button_save.config(state=tk.NORMAL)
        # remove temporary files
        os.remove(tmp_filename_in)
        os.remove(tmp_filename_out)
    
    def on_click_input_image(self, event):
        self.get_coordinates(event)
        # if coordinates are valid
        if self.x>=self.block_size and self.y>=self.block_size:
            self.overlay_block(self.prev_rect,"blue",4)
            self.root.update()
            # extract low-resolution block from given coordinates
            block = self.img_input.crop((self.x-self.block_size, self.y-self.block_size, self.x, self.y))
            try:
                self.upsample_data(block)
            except:
                print("An error occurred during FRESH upsampling...")
                print(traceback.format_exc())
    
    def upsample_image(self):
        msg = "Applying FRESH on the complete image might take a long time. Are you sure you want to proceed?"
        if tkMessageBox.askyesno(self.title, msg, icon="warning"):
            try:
                # overlay rectangle as big as the image itself
                self.scroll_img_input.get_canvas().itemconfig(self.prev_rect, outline="blue", width=4)
                self.scroll_img_input.get_canvas().coords(self.prev_rect, 2, 2, self.img_input.width-2, self.img_input.height-2)
                self.root.update()
                # run upsampling on whole image
                self.upsample_data(self.img_input)
            except:
                print("An error occurred during FRESH upsampling...")
                print(traceback.format_exc())
            
    def load_image(self):
        self.root.update()
        try:
            #self.filename_in = "/Users/Matteo/Documents/MATLAB/upsampling/hires/data/cameraman256.bmp"
            self.filename_in = tkFileDialog.askopenfilename(parent=self.root, 
                initialdir="", title='Choose an image')
            if self.filename_in is not None:
                self.img_input = Image.open(self.filename_in)
                self.scroll_img_input.set_image(self.img_input)
        except:
            tkMessageBox.showerror(self.title, "Image file type not recognized. Please select valid image file.")
            print("An error occurred while opening image file...")
        
    def save_image(self):
        filename_out = tkFileDialog.asksaveasfile(mode="w", title="Save FRESH result", defaultextension=".png")
        if filename_out is not None:
            self.block_fresh.convert("RGB").save(filename_out)
        
if __name__ == "__main__":
    root = tk.Tk()
    ws = root.winfo_screenwidth()
    hs = root.winfo_screenheight()
    w = int(ws/1.5)
    h = int(hs/1.5)
    x = int((ws-w)/2)
    y = int((hs-h)/4)
    root.geometry('{}x{}+{}+{}'.format(w,h,x,y))
    FreshGui(root).load_image()
    root.mainloop()
    