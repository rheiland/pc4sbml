from ipywidgets import Output
from IPython.display import display, HTML

class CodeTab(object):

    def __init__(self):
        # self.tab = Output(layout={'height': '600px'})
        self.tab = Output(layout={'height': 'auto'})
        # self.tab.append_display_data(HTML(filename='doc/code.html'))
        custom_code_file = open("src/custom_modules/custom.cpp", "r")
        code_str = custom_code_file.read()
        # print(code_str)
        self.tab.append_display_data(HTML('<textarea readonly rows="20" cols="110">' + code_str + '</textarea>'))
        
