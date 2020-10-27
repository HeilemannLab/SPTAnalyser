"""
@author: Johanna Rahm
Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.

Create microscope.txt file for ermine analysis with camera integration time, pixel size and localization error.
"""


class Microscope():
    def __init__(self, dt, pixel_size, error, save_dir, ym_to_nm=True):
        self.dt = dt
        self.pixel_size = pixel_size
        self.error = error * 1000 if ym_to_nm else error  # SPTAnalyser unit in ym -> convert it to nm
        self.save_dir = save_dir

    def write_microscope_file(self, file_path):
        file = open(file_path, "w+")
        if not file.closed:
            file.write("# SMLMS Microscope File \n")
            file.write("# pxl Size[nm] \n")
            file.write("# integration Time [s] \n")
            file.write("# localization precision [nm] \n")
            file.write("%.6e \n" %(float(self.pixel_size)))
            file.write("%.6e \n" %(float(self.dt)))
            file.write("%.6e \n" %(float(self.error)))

    def save_hmm_microscope(self):
        out_file_name = self.save_dir + "\\" + "microscope.txt"
        self.write_microscope_file(out_file_name)
