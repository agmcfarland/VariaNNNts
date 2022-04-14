
import subprocess
import os 
from os.path import join as pjoin
import re
from Bio.Seq import Seq
from Bio import SeqIO, SeqFeature
import glob
import shutil
import pandas as pd


pd.set_option('display.max_columns', None)


workdir = '/home/agmcfarland/flu_project/test_variant_genomes'
outdir = pjoin(workdir,'synthetic_reads')
os.makedirs(outdir,exist_ok=True)

os.chdir(workdir)


## Deletions need to be at least a start of n and stop of n+1
## Insertions can be a start of n and stop of n
## Substitutions can be a start of n and stop of n

df_master = pd.read_csv('insilico_variants.csv')


class VariantMaker:
	'''
	Make reproducible variants from a given sequence
	'''

	@classmethod
	def var_Sequence(cls, referenceSeq, df_Var):
		'''
		Takes as input a reference sequence and a dataframe of variant types and starting and stopping positions.
		Returns a sequence containign the original refernce sequence except for variants at the locations specified by df_Var
		'''
		bases = ['A','C','T','G'] # list of bases allowed for subsitutions
 
		var_positions = df_Var['STA'].tolist() # list of starting positions for all variants

		var_seq = '' # variable to store growing variant sequence 

		del_range = range(-10,-1) # intitalize a deletion range that will never be found

		base_position = 0

		for i in referenceSeq:

			base_position += 1
			
			if base_position in var_positions:
				
				## Deletions
				if df_Var[df_Var['STA']==base_position]['TYP'].item()=='DEL':
					del_range = range(df_Var[df_Var['STA']==base_position]['STA'].item(), df_Var[df_Var['STA']==base_position]['STO'].item())

				## Insertions 
				if df_Var[df_Var['STA']==base_position]['TYP'].item()=='INS':
					inserted_sequence = df_Var[df_Var['STA']==base_position]['MOD'].item()
					var_seq += i + inserted_sequence

				## Subsitutions
				if df_Var[df_Var['STA']==base_position]['TYP'].item()=='SUB':
					var_seq += [b for b in bases if b != i][0]

			if base_position not in del_range:
				
				var_seq += i

		return var_seq



df_grouped = df_master.groupby(['reference','newfile'])

print(len(df_grouped.groups))

for name, group in df_grouped:

	# reference files and filepaths don't have a suffix so that filetypes can be changed depending on output needed.
	reference_file = name[0].replace('.fasta','') # no .fasta suffix. Done for consistency
	variant_file = name[1]

	reference_filepath = pjoin(os.getcwd(),reference_file)
	variant_filepath = pjoin(outdir,variant_file)

	# the segments that will have variants introduced
	segments_to_modify = group['segment_id'].unique().tolist()

	with open(variant_filepath+'.fasta', 'w') as outfile:

		for record in SeqIO.parse(reference_filepath+'.fasta', 'fasta'):
			_segment_id = int(record.id[-1:])


			if _segment_id in segments_to_modify:

				df_segment_variants = group[group['segment_id']==_segment_id][['TYP','STA','STO','MOD']]

				writeseq = VariantMaker.var_Sequence(referenceSeq = str(record.seq), df_Var = df_segment_variants)
			else:
				writeseq = str(record.seq)

			outfile.write('>{}\n{}\n'.format(record.id, writeseq))

	cmd = [
		'randomreads.sh',
		'reads=500000',
		# 'coverage=1000',
		'prefix={}'.format(variant_file),
		'seed=5','paired=t','addpairnum=f','illuminanames=t','adderrors=f','usequality=false',
		'ref={}.fasta'.format(variant_filepath),'out1={}_R1_001.fastq'.format(variant_filepath),'out2={}_R2_001.fastq'.format(variant_filepath),
		'minlength=150','maxlength=150',
		'maxsnps=0',
		'maxinss=0',
		'maxdels=0',
		'maxsubs=0',
		'q=40',
		'superflat=t',
		'overwrite=t',
		'qv=0']

	subprocess.run(cmd)


os.system('rm -R {}/*.gz'.format(outdir))
os.system('parallel --gnu gzip  ::: {}/*.fastq'.format(outdir))















