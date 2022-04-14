import sys
import os 
from os.path import join as pjoin
from datetime import datetime
import logging

class RunLogger:
	'''
	Handles logging: Making a unique logfile by name. Adding messages to it. Removing logfiles.
	'''
	def __init__(self, directory, filename):
		'''
		directory is where the logfile will be stored
		filename is the name of the file. no suffix is supplied
		'''
		self.directory = directory
		self.filename = filename
		self.filepath = pjoin(self.directory,self.filename)
		self.formatter = logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(message)s',
			datefmt='%Y-%m-%d %H:%M:%S')

	def initialize_FileHandler(self):
		'''
		Initialize a logfile with time stamps and in append mode.
		'''
		handler = logging.FileHandler('{}.txt'.format(self.filepath), mode='a')
		handler.setFormatter(self.formatter)
		self.logger = logging.getLogger(self.filename)
		self.logger.setLevel(logging.DEBUG)
		self.logger.addHandler(handler)

	def add_StreamHandler(self):
		'''
		Adds the stream handler to the logger object
		'''
		handler = logging.StreamHandler()
		handler.setFormatter(self.formatter)
		self.logger.addHandler(handler)