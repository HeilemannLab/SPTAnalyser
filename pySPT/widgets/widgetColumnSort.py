# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 15:52:19 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""


class WidgetColumnSort():
    def __init__(self, file_name, file_type, significant_words):
        #self.file = []
        self.file_name = file_name
        if file_type == "rapidSTORM":
            self.identifier_after = "syntax"  # identifier after target word
        elif file_type == "swift":
            self.identifier_after = "unit"
        self.significant_words = significant_words
        #self.significant_words = ['"track_id"', '"mjd_n"', '"mjd"']
        self.identifier_before = "identifier"  # identifier before target word
        self.header = ""
        self.number_columns = 0
        self.sub_headers = [] # index in list = index of column in file
        self.target_words= []  # raw target words in double quotes
        self.column_order = {}    
        
    def load_file(self):
        file = open(self.file_name)
        self.header = file.readline()  # get the header as first line

    def testing_header(self):
        """
        Test different headers without loading.
        """
        #self.header =  '<localizations insequence="true" repetitions="variable"><field identifier="mjd" unit=""/><field identifier="track_id" unit=""/><field identifier="mjd_n" unit=""/></localizations>'
        self.header = '<localizations insequence="true" repetitions="variable"><field identifier="Position-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position in sample space in X" unit="nanometer" min="0 m" max="4.029e-005 m" /><field identifier="Position-0-0-uncertainty" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position uncertainty in sample space in X" unit="nanometer" /><field identifier="Position-1-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position in sample space in Y" unit="nanometer" min="0 m" max="4.029e-005 m" /><field identifier="Position-1-0-uncertainty" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position uncertainty in sample space in Y" unit="nanometer" /><field identifier="ImageNumber-0-0" syntax="integer" semantic="frame number" unit="frame" min="0 fr" max="999 fr" /><field identifier="Amplitude-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="emission strength" unit="A/D count" /><field identifier="PSFWidth-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="PSF FWHM in X" unit="nanometer" /><field identifier="PSFWidth-1-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="PSF FWHM in Y" unit="nanometer" /><field identifier="FitResidues-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="fit residue chi square value" unit="dimensionless" /><field identifier="LocalBackground-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="local background" unit="A/D count" /></localizations>'
                
    def sub_header(self):
        """
        Create list with sub headers in order of columns.
        """
        cut_header = self.header
        while cut_header.find(self.identifier_before) != -1:
            slice_index = cut_header.find(self.identifier_before)
            sub_header = cut_header[:slice_index+len(self.identifier_before)]
            self.sub_headers.append(sub_header)
            cut_header = cut_header[slice_index+len(self.identifier_before):]
        self.sub_headers.append(cut_header)
        self.sub_headers.pop(0)
        self.number_columns = len(self.sub_headers)
        
    def get_words(self):
        #cut_header = self.sub_headers
        for i in range(0, len(self.sub_headers)):
            sub_header = self.sub_headers[i]
            word_end_index = sub_header.find(self.identifier_after) - 1
            target_word = sub_header[1:word_end_index]
            self.target_words.append(target_word)
    
    def column_index(self):
        for target_word in self.target_words:
            for word in self.significant_words:
                if word in target_word:
                    if word not in self.column_order:
                    #append word (value) and index of sub_head (key) to dictionary
                        self.column_order[self.target_words.index(target_word)] = word
        
    def run_column_sort(self):
        self.load_file()
        self.sub_header()
        self.get_words()
        self.column_index()

        
def main():
    file_name = "F:/Marburg/single_colour_tracking/resting/160404_CS5_Cell1/cell_1_MMStack_Pos0.ome.tif.txt"
    widget_column_sort = WidgetColumnSort(file_name, "rapidSTORM", ['"Position-0-0-uncertainty"', '"Position-1-0-uncertainty"'])
    widget_column_sort.run_column_sort()
    print(widget_column_sort.column_order)
    
    
if __name__ == "__main__":
    main()
