# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 15:52:19 2019

@author: pcoffice37

Research group Heilemann
Institute for Physical and Theoretical Chemistry, Goethe University Frankfurt a.M.
"""


class ColumnSort():
# =============================================================================
#     def __init__(self, file_type, significant_words):
#         #self.file = []
#         self.file_name = ""
#         if file_type == "rapidSTORM":
#             self.word_identifier = "syntax"  # next word after target word (can be syntax, unit ...)
#         elif file_type == "swift":
#             self.word_identifier = "unit"
#         self.significant_words = significant_words
#         #self.significant_words = ['"track_id"', '"mjd_n"', '"mjd"']
#         self.identifier = "identifier"  # identifier
#         self.header = ""
#         self.number_columns = 0
#         self.sub_headers = [] # index in list = index of column in file
#         self.target_words= []  # raw target words in double quotes
#         self.column_order = {}      
# =============================================================================
    def __init__(self):
        #self.file = []
        self.file_name = ""
        self.word_identifier = "syntax"  # next word after target word (can be syntax, unit ...)
        self.significant_words = ['"Position-0-0"', '"PSFWidth-1-0"', '"Position-1-0-uncertainty"', '"Position-1-0"']
        self.identifier = "identifier"  # identifier
        self.header = ""
        self.number_columns = 0
        self.sub_headers = [] # index in list = index of column in file
        self.target_words= []  # raw target words in double quotes
        self.column_order = {}   
# =============================================================================
#     def load_file(self, file_name):
#         file = open(file_name)
#         self.header = file.readline()  # get the header as first line
#         print(self.header)
# =============================================================================
    
    def testing_header(self):
        
        #self.header =  '<localizations insequence="true" repetitions="variable"><field identifier="mjd" unit=""/><field identifier="track_id" unit=""/><field identifier="mjd_n" unit=""/></localizations>'
        self.header = '<localizations insequence="true" repetitions="variable"><field identifier="Position-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position in sample space in X" unit="nanometer" min="0 m" max="4.029e-005 m" /><field identifier="Position-0-0-uncertainty" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position uncertainty in sample space in X" unit="nanometer" /><field identifier="Position-1-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position in sample space in Y" unit="nanometer" min="0 m" max="4.029e-005 m" /><field identifier="Position-1-0-uncertainty" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position uncertainty in sample space in Y" unit="nanometer" /><field identifier="ImageNumber-0-0" syntax="integer" semantic="frame number" unit="frame" min="0 fr" max="999 fr" /><field identifier="Amplitude-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="emission strength" unit="A/D count" /><field identifier="PSFWidth-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="PSF FWHM in X" unit="nanometer" /><field identifier="PSFWidth-1-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="PSF FWHM in Y" unit="nanometer" /><field identifier="FitResidues-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="fit residue chi square value" unit="dimensionless" /><field identifier="LocalBackground-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="local background" unit="A/D count" /></localizations>'

        print(self.header)
        
# =============================================================================
#     def sub_header(self):
#         pass
#         cut_header = self.header
#         while cut_header.find(self.identifier) != -1:
#             slice_index = cut_header.find(self.identifier)
#             print(slice_index)
#             sub_header = cut_header[:slice_index+len(self.identifier)]
#             self.sub_headers.append(sub_header)
#             cut_header = cut_header[slice_index+len(self.identifier):]
#             #print(sub_header)
#             #print(self.sub_headers)
#             #print(cut_header)
#         self.sub_headers.append(cut_header)
#         self.sub_headers.pop(0)
#         print(self.sub_headers)
#         self.number_columns = len(self.sub_headers)
# =============================================================================
        
    def sub_header(self):
        cut_header = self.header
        while cut_header.find(self.identifier) != -1:
            slice_index = cut_header.find(self.identifier)
            print(slice_index)
            sub_header = cut_header[:slice_index+len(self.identifier)]
            self.sub_headers.append(sub_header)
            cut_header = cut_header[slice_index+len(self.identifier):]
            #print(sub_header)
            #print(self.sub_headers)
            #print(cut_header)
        self.sub_headers.append(cut_header)
        self.sub_headers.pop(0)
        print(self.sub_headers)
        self.number_columns = len(self.sub_headers)
        #for i in self.sub_headers:
        #    i.strip()
        #print(self.sub_headers)
        
    def get_words(self):
        #cut_header = self.sub_headers
        for i in range(0, len(self.sub_headers)):
            sub_header = self.sub_headers[i]
            word_end_index = sub_header.find(self.word_identifier) - 1
            target_word = sub_header[1:word_end_index]
            self.target_words.append(target_word)
        print(self.target_words)
        
# =============================================================================
#     def column_index(self):
#         for sub_header in self.sub_headers:
#             #word_found = False
#             for word in self.significant_words:
# 
#                 #while not word_found:
#                 if word in sub_header:
#                     if word not in self.column_order:
#                     #append word (value) and index of sub_head (key) to dictionary
#                         self.column_order[self.sub_headers.index(sub_header)] = word
#                     #word_found = True
#         print(self.column_order)
# =============================================================================
            
    
    def column_index(self):
        for target_word in self.target_words:
            #word_found = False
            for word in self.significant_words:
                print("Target word, word", target_word, word)

                #while not word_found:
                if word in target_word:
                    if word not in self.column_order:
                    #append word (value) and index of sub_head (key) to dictionary
                        self.column_order[self.target_words.index(target_word)] = word
                    #word_found = True
        print("Dict", self.column_order)
    
    def get_key(self):
        pass
        
        #print(self.column_order.keys()[self.column_order.index("track_id")])
        
        print(type(list(self.column_order.keys())[list(self.column_order.values()).index('"track_id"')]))

      
        
def main():
    column_sort = ColumnSort()
    #file_name = "x"
    #file_name = "xidentifiertrack_ididentifiermjd_nfkkaidentifiermjd"
    #sorting_columns.header = file_name
    
    #file_name = "F:\\Marburg\\single_colour_tracking\\resting\\160404_CS5_Cell1\\cell_1_MMStack_Pos0.ome.tif.tracked.seg.txt"

    #sorting_columns.load_file(file_name)
    column_sort.testing_header()
    column_sort.sub_header()
    column_sort.get_words()
    column_sort.column_index()
    #sorting_columns.get_key()
    
if __name__ == "__main__":
    main()
