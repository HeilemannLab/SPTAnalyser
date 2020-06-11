# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 15:52:19 2019

@author: Johanna Rahm

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""


class WidgetColumnSort():
    def __init__(self, file_name, file_type, significant_words):
        self.file_name = file_name
        if file_type == "rapidSTORM":
            self.software = file_type
            self.identifier_after = "syntax"
        elif file_type == "ThunderSTORM":
            self.software = file_type
            self.identifier_after = ","
        self.significant_words = significant_words
        self.identifier_before = "identifier"  # identifier before target word
        self.header = ""
        self.number_columns = 0
        self.sub_headers = [] # index in list = index of column in file
        self.target_words= []  # raw target words in double quotes
        self.column_order = {}  # {0: '"track_id"', 4: '"mjd"', 6: '"mjd_n"'}
        self.correct_header = False
        
    def check_header(self):
        """
        Open header & check if significant words are in header, if they
        appear count will go up 1, if word appears multible times it does not
        matter, count will still go up only once.
        """
        file = open(self.file_name)
        self.header = file.readline()  # get the header as first line
        count = 0
        for significant_word in self.significant_words:  # iterate through list of target columns
            for i in significant_word:  # iterate through tuples (tuple = different possible names for same column)
                if i in self.header:
                    count += 1
        if count == len(self.significant_words):
            self.correct_header = True

    def testing_header(self):
        """
        Test different headers without loading.
        """
        #self.header =  '<localizations insequence="true" repetitions="variable"><field identifier="mjd" unit=""/><field identifier="track_id" unit=""/><field identifier="mjd_n" unit=""/></localizations>'
        self.header = '<localizations insequence="true" repetitions="variable"><field identifier="Position-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position in sample space in X" unit="nanometer" min="0 m" max="4.029e-005 m" /><field identifier="Position-0-0-uncertainty" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position uncertainty in sample space in X" unit="nanometer" /><field identifier="Position-1-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position in sample space in Y" unit="nanometer" min="0 m" max="4.029e-005 m" /><field identifier="Position-1-0-uncertainty" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position uncertainty in sample space in Y" unit="nanometer" /><field identifier="ImageNumber-0-0" syntax="integer" semantic="frame number" unit="frame" min="0 fr" max="999 fr" /><field identifier="Amplitude-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="emission strength" unit="A/D count" /><field identifier="PSFWidth-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="PSF FWHM in X" unit="nanometer" /><field identifier="PSFWidth-1-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="PSF FWHM in Y" unit="nanometer" /><field identifier="FitResidues-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="fit residue chi square value" unit="dimensionless" /><field identifier="LocalBackground-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="local background" unit="A/D count" /></localizations>'
                
    def ts_sub_headers(self):
        self.sub_headers = self.header.split(",")
        self.sub_headers[-1] = self.sub_headers[-1][:-1]  # get rid of the new line character for the last sub_header

    def ts_create_column_order(self):
        """
        thunderSTORM: Add sub header index and value to dictionary.
        """
        for i in self.sub_headers:
            for significant_word in self.significant_words:
                if i in significant_word:
                    self.column_order.update({self.sub_headers.index(i):i})
        
    def rs_sub_headers(self):
        """
        rapidSTORM: Create list with sub headers in order of columns.
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
        
    def rs_get_words(self):
        """
        rapidSTORM
        """
        #cut_header = self.sub_headers
        for i in range(0, len(self.sub_headers)):
            sub_header = self.sub_headers[i]
            word_end_index = sub_header.find(self.identifier_after) - 1
            target_word = sub_header[1:word_end_index]
            self.target_words.append(target_word)
    
    def rs_column_index(self):
        """
        rapidSTORM
        """
        for target_word in self.target_words:
            for words in self.significant_words:
                for word in words:
                    if word in target_word:
                        if word not in self.column_order:
                        #append word (value) and index of sub_head (key) to dictionary
                            self.column_order[self.target_words.index(target_word)] = word
        
    def run_column_sort(self):
        if self.software == "ThunderSTORM":
            self.ts_sub_headers()
            self.ts_create_column_order()
        elif self.software == "rapidSTORM":
            self.rs_sub_headers()
            self.rs_get_words()
            self.rs_column_index()
        
        
def main():
    file_name = "F:/Marburg/single_colour_tracking/resting/160404_CS5_Cell1/cell_1_MMStack_Pos0.ome.tif.txt"
    widget_column_sort = WidgetColumnSort(file_name, "rapidSTORM", ['"Position-0-0-uncertainty"', '"Position-1-0-uncertainty"'])
    widget_column_sort.run_column_sort()
    print(widget_column_sort.column_order)
    
    
if __name__ == "__main__":
    main()
