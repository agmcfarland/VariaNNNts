import os
import sys
import subprocess
from .RunLogger import RunLogger

def call_Command(cmd, logger_, shell_=False):
	'''
	Runs a shell command using subprocess.run. Recods output to a logger object that has already been created.
	Default is to use subprocess.run. If shell is necessary then subprocess.Popen is used.
	'''
	try:
		if shell_ == False:
			capture = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True).stdout
			logger_.logger.info(capture)
		else:
			capture = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			logger_.logger.info(capture.communicate()[0].decode('utf-8'))
	except:
		logger_.logger.exception('\nError:', exc_info=True)
