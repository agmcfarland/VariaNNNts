import pandas as pd
import numpy as np

class VariantMaker:
	'''
	Make reproducible variants from a given sequence
	'''
	@classmethod
	def var_Sequence(cls, referenceSeq, df_Var):
		'''
		Takes as input a reference sequence and a pandas dataframe of variant types and starting and stopping positions.
		Returns a sequence containing the original refernce sequence except for variants at the locations specified by df_Var
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