import os
import xml.sax

path = "./TC_ILP/"
outfile = open(path+"Sol_summary.txt","w")

class solHandler(xml.sax.ContentHandler):
	def __init__(self):
		self.content = ""
		self.tag = ""

	def startElement(self,tag,attributes):
		self.tag = tag
		if self.tag == "header":
			outfile.write(attributes["problemName"]+"\n")
		if self.tag == "variable":
			if attributes["name"] == "M1":
				outfile.write("("+attributes["value"]+" ")# MUL:
			elif attributes["name"] == "M2":
				outfile.write(attributes["value"]+")\n")# ALU:

	def endElement(self,tag):
		pass

	def characters(self,content):
		self.content = content

for num in ["1.0","1.5","2.0"]:		
	for file in os.listdir(path):
		if (file[-4:] == ".sol" and num in file):
			print("Start parsing "+file+"...")
			infile = open(path+file,"r")
			parser = xml.sax.make_parser()
			handler = solHandler()
			parser.setContentHandler(handler)
			parser.parse(path+file)
			print("Finish parsing "+file+".")

outfile.close()