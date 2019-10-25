from ipywidgets import Output
from IPython.display import display, HTML, FileLink
from display_xml import XML

class DefaultXMLTab(object):

    def __init__(self):
        # self.tab = Output(layout={'height': '600px'})
        self.tab = Output(layout={'height': 'auto'})
        xml_file = open("data/PhysiCell_settings.xml", "r")
        xml_str = xml_file.read()
        # print(xml_str)
        self.tab.append_display_data(HTML('<textarea readonly rows="20" cols="110">' + xml_str + '</textarea>'))
        # self.tab.append_display_data(HTML('<textarea rows="20" cols="40" style="border:none;">' + '<fred>42</fred>\n' + '</<textarea rows="20" cols="40" style="border:none;">>\n'))
        # self.tab.append_display_data(HTML('<textarea style="border:none;">' + '<fred>42</fred>\n' + '</textarea>\n'))
        # self.tab.append_display_data(XML('<pre>' + xml_str + '</pre>'))
        # self.tab.append_display_data(XML(xml_str))
        # with self.tab:
        #     display(HTML(filename='data/PhysiCell_settings.xml'))
        
