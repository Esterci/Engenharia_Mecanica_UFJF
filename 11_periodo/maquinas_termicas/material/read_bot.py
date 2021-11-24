from io import StringIO
from os import read
from pdfminer.converter import TextConverter
from pdfminer.layout import LAParams
from pdfminer.pdfdocument import PDFDocument
from pdfminer.pdfinterp import PDFResourceManager, PDFPageInterpreter
from pdfminer.pdfpage import PDFPage
from pdfminer.pdfparser import PDFParser
import PyPDF2
import argparse


def convert_pdf_to_string(file_path):

	output_string = StringIO()
	with open(file_path, 'rb') as in_file:
	    parser = PDFParser(in_file)
	    doc = PDFDocument(parser)
	    rsrcmgr = PDFResourceManager()
	    device = TextConverter(rsrcmgr, output_string, laparams=LAParams())
	    interpreter = PDFPageInterpreter(rsrcmgr, device)
	    for page in PDFPage.create_pages(doc):
	        interpreter.process_page(page)

	return(output_string.getvalue())


parser = argparse.ArgumentParser(description = '', add_help = False)
parser = argparse.ArgumentParser()


parser.add_argument('-f','--file', action='store',
        dest='file', required = True,
            help = "The PDF file to be converted.")


args = parser.parse_args()

file = args.file

file_name = file.split('.')[0]

reader = PyPDF2.PdfFileReader(file)

print(reader.documentInfo)

num_of_pages = reader.numPages

print('Number of pages: ' + str(num_of_pages))

text = convert_pdf_to_string(file)

text = text.replace('\x0c',' ')
text = text.replace('%','\%')
formated_text = text.split('\n')

for line in formated_text:

    # Open the file in append & read mode ('a+')
    with open(file_name+"_conv.txt", "a+") as file_object:
        # Move read cursor to the start of file.
        file_object.seek(0)
        # If file is not empty then append '\n'
        data = file_object.read(100)
        if len(data) > 0 :
            file_object.write("\n")
        # Append text at the end of file
        file_object.write(line)

print('Convertion complete!')