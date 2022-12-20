from io import StringIO
from os import read
import argparse


class ProgBar:

    def __init__(self, n_elements,int_str):
        
        import sys

        self.n_elements = n_elements
        self.progress = 0

        print(int_str)

        # initiallizing progress bar

        info = '{:.2f}% - {:d} of {:d}'.format(0,0,n_elements)

        formated_bar = ' '*int(50)

        sys.stdout.write("\r")

        sys.stdout.write('[%s] %s' % (formated_bar,info))

        sys.stdout.flush()

    def update(self,prog_info=None):
        
        import sys

        if prog_info == None:

            self.progress += 1

            percent = (self.progress)/self.n_elements * 100 / 2

            info = '{:.2f}% - {:d} of {:d}'.format(percent*2,self.progress,self.n_elements)

            formated_bar = '-'* int (percent) + ' '*int(50-percent)

            sys.stdout.write("\r")

            sys.stdout.write('[%s] %s' % (formated_bar,info))

            sys.stdout.flush()


        else:

            self.progress += 1

            percent = (self.progress)/self.n_elements * 100 / 2

            info = '{:.2f}% - {:d} of {:d} '.format(percent*2,self.progress,self.n_elements) + prog_info

            formated_bar = '-'* int (percent) + ' '*int(50-percent)

            sys.stdout.write("\r")

            sys.stdout.write('[%s] %s' % (formated_bar,info))

            sys.stdout.flush()


class PDF_miner:

    def __init__(self,folder):

        import PyPDF2

        self.folder = folder

        reader = PyPDF2.PdfFileReader(folder + '/slide.pdf')

        self.reader = PyPDF2.PdfFileReader(folder + '/slide.pdf')

        self.pdf_n_pages = reader.numPages

        self.pic_ref_list = []

        self.discart_list = []


    def convert_pdf_to_string(self):

        from io import StringIO
        from pdfminer.converter import TextConverter
        from pdfminer.layout import LAParams
        from pdfminer.pdfdocument import PDFDocument
        from pdfminer.pdfinterp import PDFResourceManager, PDFPageInterpreter
        from pdfminer.pdfpage import PDFPage
        from pdfminer.pdfparser import PDFParser


        output_string = StringIO()
        with open('transition.pdf', 'rb') as in_file:
            parser = PDFParser(in_file)
            doc = PDFDocument(parser)
            rsrcmgr = PDFResourceManager()
            device = TextConverter(rsrcmgr, output_string, laparams=LAParams())
            interpreter = PDFPageInterpreter(rsrcmgr, device)
            for page in PDFPage.create_pages(doc):
                interpreter.process_page(page)

        return(output_string.getvalue())


    def add_line(self,line,out):

        # Open the file in append & read mode ('a+')

        with open(self.folder+out, "a+") as file_object:

            # Move read cursor to the start of file.
            file_object.seek(0)

            # If file is not empty then append '\n'
            data = file_object.read(100)

            if len(data) > 0 :
                file_object.write("\n")
                
            # Append text at the end of file
            file_object.write(line)

    
    def get_text(self,pg):
        
        import PyPDF2

        writer = PyPDF2.PdfFileWriter()

        output_filename = 'transition.pdf'

        page = self.reader.getPage(pg)

        writer.addPage(page)

        with open(output_filename, 'wb') as output:
            writer.write(output)

        del writer

        text = self.convert_pdf_to_string()

        text = text.replace('\x0c',' ')
        text = text.replace('%','\%')
        text = text.split('\n')
        # text = text[5:]

        return text


    def create_fig_form(self,type="fig"):

        for i in range(1,self.pdf_n_pages+1):

            if type == "fig": name = "/fig_form.txt"

            else: name = "/eq_form.txt"

            self.add_line("{}:0".format(i),name)


    def extract_figures(self):

        import fitz
        import glob
        
        print('Commencing Figure Extraction')
        
        
        # opening pdf file

        self.doc = fitz.open(self.folder + "/slide.pdf")

        figures_list = glob.glob(self.folder + '/figures/*')

        for pg in range(self.pdf_n_pages):

            for img in self.doc.get_page_images(pg):
                
                # get figure refferance
                
                xref = img[0]
                
                # build figure bitmap

                pix = fitz.Pixmap(self.doc, xref)

                # defines the name of the picture

                file_name = self.folder + ("/figures/%s.png" % (xref))

                # if the figure does not exists save it

                if not file_name in figures_list:

                    if pix.n < 5:       # this is GRAY or RGB
                        pix.save(file_name)
                    else:               # CMYK: convert to RGB first
                        pix1 = fitz.Pixmap(fitz.csRGB, pix)
                        pix1.save(file_name)
                        pix1 = None
                    pix = None

        print('Figure Extraction done!')

        print('='*22)

        self.create_fig_form()

        print("\n")

        input("Fill th figures location form and press ENTER")

        with open(self.folder + "/fig_form.txt", "r+") as file_object:

            fig_loc = file_object.read()

        fig_loc = fig_loc.replace(' ','')
        fig_loc = fig_loc.split('\n')

        self.fig_loc = fig_loc


    def extract_equations(self):

        import numpy as np
        import shutil
        import glob

        print("\nPress ENTER to confirm that the equation prints are in the correct folder")

        input("Press ENTER to confirm that the equation prints are in the correct folder")

        figures_list = glob.glob(self.folder + '/equations/*')

        figures_list = np.sort(figures_list)

        for i,fig in enumerate(figures_list):

            file_name = self.folder + '/equations/' + str(i+1) + '.png'

            shutil.move(fig, file_name)

        self.create_fig_form("eq")

        print("\nFill the equations location form and press ENTER")

        input("Fill the equations location form and press ENTER")

        with open(self.folder + "/eq_form.txt", "r+") as file_object:

            fig_loc = file_object.read()

        fig_loc = fig_loc.replace(' ','')
        fig_loc = fig_loc.split('\n')

        self.eq_loc = fig_loc


    def get_figures (self,pg):

        pg_figs = self.fig_loc[pg]

        pg_figs = pg_figs.split(':')[1:]

        pg_figs = pg_figs[0].split(',')

        fig_call_list = []

        for fig in pg_figs:

            if not fig == '0':

                figure_call = [
                    '\\begin{figure}[h!]',
                    '\centering',
                    ('\includegraphics[width=0.8\\textwidth]{Pictures/' + fig + '.png}'),
                    '\caption*{}',
                    '\label{fig:my_label}',
                    '\end{figure}',
                    '',
                    ''
                ]

                for i in range(len(figure_call)): fig_call_list.append(figure_call[i])


        return fig_call_list


    def get_equations(self,pg):

        pg_eq = self.eq_loc[pg]

        pg_eq = pg_eq.split(':')[1:]

        pg_eq = pg_eq[0].split(',')

        fig_call_list = []

        for fig in pg_eq:

            if not fig == '0':

                figure_call = [
                    '\\begin{figure}[h!]',
                    '\centering',
                    ('\includegraphics[width=0.6\\textwidth]{Equations/' + fig + '.png}'),
                    '\caption*{}',
                    '\label{fig:my_label}',
                    '\end{figure}',
                    '',
                    ''
                ]

                for i in range(len(figure_call)): fig_call_list.append(figure_call[i])


        return fig_call_list


    def convert(self):

        self.pages_struct = {}

        self.extract_figures()

        self.extract_equations()

        # Initializing progress bar

        bar = ProgBar(self.pdf_n_pages,'\nConverting PDF to LaTex...')

        for pg in range(self.pdf_n_pages):

            self.pages_struct[pg] = {}

            self.pages_struct[pg]['text'] = self.get_text(pg)

            self.pages_struct[pg]['equations'] = self.get_equations(pg)

            self.pages_struct[pg]['figures'] = self.get_figures(pg)

            bar.update()


    def save(self):

        ### Saving results in a txt file 

        # initializing progress bar

        bar = ProgBar(len(self.pages_struct),"\n\nSalving conversion...")

        for pg in self.pages_struct:
            
            # appending text

            for txt in self.pages_struct[pg]['text']:
                self.add_line(txt,"/conv.txt")

            # appending equations

            for fig in self.pages_struct[pg]['equations']:
                self.add_line(fig,"/conv.txt")

            # appending figures

            for fig in self.pages_struct[pg]['figures']:
                self.add_line(fig,"/conv.txt")

            # adding new page

            self.add_line('\\newpage',"/conv.txt")

            bar.update()
  
parser = argparse.ArgumentParser(description = '', add_help = False)
parser = argparse.ArgumentParser()


parser.add_argument('-f','--folder', action='store',
        dest='folder', required = True,
            help = "The folder to perform the conversion.")


args = parser.parse_args()

folder = args.folder

# initializing class

converter = PDF_miner(folder)

# converting information

converter.convert()

# saving results


converter.save()