#!/usr/bin/env python
import sys
import subprocess
import os 
from os.path import join as pjoin
from Bio.Seq import Seq
from Bio import SeqIO
import shutil
import pandas as pd
import psutil
import argparse

global script_path
script_path = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0]
sys.path.append(script_path)

from libs.Misc import call_Command
from libs.RunLogger import RunLogger
from libs.VariantMaker import VariantMaker


def generate_VariantSequences(args):
	'''
	Creates an assembly and synthetic reads contaning specified variants
	'''
	if args.overwrite == True:
		try:
			shutil.rmtree(args.output_directory)
		except:
			pass

	os.makedirs(args.output_directory, exist_ok=True)

	# Start logging. Record start, memory available, threads available, and arguments.
	run_logger = RunLogger(directory=args.output_directory,filename='runLog')
	run_logger.initialize_FileHandler()
	run_logger.add_StreamHandler()
	run_logger.logger.info('Starting VariaNNNts...\n')
	run_logger.logger.info('Total memory available: {}Gb\n'.format(round(psutil.virtual_memory()[0]/(1024.0**3)),2))
	run_logger.logger.info('Threads available: {}\n'.format(os.cpu_count()))
	run_logger.logger.info('Program arguments:\n')
	arguments_list = vars(args)
	for k,v in arguments_list.items():
		run_logger.logger.info('{}: {}\n'.format(k,v))


	# read in variant table file. if it doesn't exist error out
	run_logger.logger.info('Reading in variant file...\n')
	try:
		df_master = pd.read_csv(args.variant_table)
		df_grouped = df_master.groupby(['reference','newfile'])
		run_logger.logger.info('Making {} new genomes and reads with variants.\n'.format(len(df_grouped.groups)))

	except:
		sys.exit(run_logger.logger.exception('No variant file could be read.\n'))


	run_logger.logger.info('Making variant reads...\n')
	# Make variant sequences	
	for name, group in df_grouped:

		# reference files and filepaths don't have a suffix so that filetypes can be changed depending on output needed.
		reference_file = name[0].replace('.fasta','') # no .fasta suffix. Done for consistency
		variant_file = name[1]
		reference_filepath = pjoin(args.genome_directory,reference_file)
		variant_filepath = pjoin(args.output_directory,variant_file)

		# the segments that will have variants introduced
		segments_to_modify = group['segment_id'].unique().tolist()

		run_logger.logger.info('{}: making genome\n'.format(variant_file))
		with open(variant_filepath+'.fasta', 'w') as outfile:

			for record in SeqIO.parse(reference_filepath+'.fasta', 'fasta'):
				_segment_id = int(record.id[-1:])


				if _segment_id in segments_to_modify:

					df_segment_variants = group[group['segment_id']==_segment_id][['TYP','STA','STO','MOD']]

					writeseq = VariantMaker.var_Sequence(referenceSeq = str(record.seq), df_Var = df_segment_variants)
				else:
					writeseq = str(record.seq)

				outfile.write('>{}\n{}\n'.format(record.id, writeseq))

		run_logger.logger.info('{}: making reads\n'.format(variant_file))
		call_Command(cmd=
			['randomreads.sh',
			'reads=500000',
			# 'coverage=1000',
			'prefix={}'.format(variant_file),
			'seed={}'.format(args.seed_bbt),
			'paired=t',
			'addpairnum=f',
			'illuminanames=t',
			'adderrors={}'.format(args.adderrors_bbt),
			'ref={}.fasta'.format(variant_filepath),'out1={}_R1_001.fastq'.format(variant_filepath),'out2={}_R2_001.fastq'.format(variant_filepath),
			'minlength=150','maxlength=150',
			'maxsnps=0',
			'maxinss=0',
			'maxdels=0',
			'maxsubs=0',
			'q={}'.format(args.q_bbt),
			'superflat=t',
			'overwrite=t',
			'qv={}'.format(args.qv_bbt)],
			logger_=run_logger)


	run_logger.logger.info('Gzipping all fastq files...\n')

	call_Command(cmd=
		'rm -R {}/*.gz'.format(args.output_directory),
		logger_=run_logger,
		shell_=True)

	call_Command(cmd=
		'parallel --gnu gzip  ::: {}/*.fastq'.format(args.output_directory),
		logger_=run_logger,
		shell_=True)

	run_logger.logger.info('Copying variant table to output directory...\n')
	shutil.copyfile(args.variant_table, pjoin(args.output_directory,os.path.basename(args.variant_table)))

	run_logger.logger.info('Finished running VariaNNNTs.\n')
	run_logger.logger.info('Output files stored in {}\n'.format(args.output_directory))


def main(args=None):
	'''
	'''
	if args is None:
		args = sys.argv[1:]

	parser = argparse.ArgumentParser(prog = 'VariaNNNts')

	parser.add_argument('--output_directory', type = str, default = pjoin(os.getcwd(),'VariaNNNts_output'), help = 'directory to output files. [./VariaNNNts_output]', metavar = '')
	parser.add_argument('--variant_table', type = str, default = pjoin(os.getcwd(),'variant_table.csv'), help = 'CSV table with desired variants. [./variant_table.csv]', metavar = '')
	parser.add_argument('--genome_directory', type = str, default = os.getcwd(), help = 'directory containing genomes to generate variants from. [./]', metavar = '')
	parser.add_argument('--overwrite', action = 'store_true', default=False, help = 'flag to overwrite all files in the output_directory. [False]')
	parser.add_argument('--read_count_bbt', type = int, default= 500000, help='number of synthetic read pairs to generate. [500000]', metavar = '')
	parser.add_argument('--seed_bbt', type = int, default= 5, help='random number generator seed to make synthetic reads. [5]', metavar = '')
	parser.add_argument('--q_bbt', type = int, default = 40, help='average base quality value. use -1 for a random seed. [40]', metavar = '')
	parser.add_argument('--qv_bbt', type = int, default = 0, help='standard deviation of variation around q_bbt [0]', metavar = '')
	parser.add_argument('--adderrors_bbt', type = str, default='f', choices=['t','f'], help='add substitution errors based on quality values [f]', metavar = '')

	args = parser.parse_args()

	generate_VariantSequences(args = args)


if __name__ == '__main__':
	sys.exit(main())















