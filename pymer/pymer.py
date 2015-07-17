from Bio.SeqUtils import MeltingTemp as mt

import argparse, rnafold, math


def analyse(sequence,
						Tm=62.0,
						Na=50.,
						K=0.,
						Tris=0.,
						Mg=2.0,
						dNTPs=0.2,
						oligo=500.,
						min_length=10,
						max_length=50):
	"""Return a list of secondary structure temparatures for each position
			[(3'base, primer_len, worst_deltaG, Tm), ...]
	"""
	ret = []
	for i, base in enumerate(sequence):
		#find the appropriate length of the primer
		ok=False
		the_Tm = Tm - 100.
		for primer_length in range(min_length, min(i, max_length)):
			pseq = sequence[i-primer_length:i]
			last_Tm = mt.Tm_NN(pseq,
												 Na=Na,
												 K=K,
												 Tris=Tris,
												 Mg=Mg,
												 dNTPs=dNTPs,
												 dnac1=oligo,
												 dnac2=0.)
			if last_Tm >= Tm:
				#the last temparature was closer
				if abs(last_Tm - Tm) < abs(the_Tm - Tm):
					the_Tm = last_Tm
					primer_length = primer_length -1
					pseq = sequence[i-primer_length:i]
				ok=True
				break
		if not ok:
			ret.append((base, 0, float('-inf'), None))
			continue

		dG = rnafold.fold(pseq)

		ret.append((base, primer_length, dG, the_Tm))

	return ret

rc = {'A':'T','T':'A','G':'C','C':'G'}
def reverse_complement(seq):
	return ''.join(rc[x] for x in reversed(seq.upper()))



def is_valid_sequence(parser, x):
	x = str(x).upper()
	invalid = []
	for char in x:
		if char not in ['A','T','C','G',]:
			invalid.append(char)
	if invalid:
		parser.error(('Invalid character{}: \'{}\', '+
								  'only \'A\', \'T\', \'G\', \'C\'')
								  .format('s' if len(invalid) > 1  else'',
								 				  '\', \''.join(invalid)))
	else:
		return x

def print_output(sequence, fwd, rev):
	rsequence = reverse_complement(sequence)
	max_flength = max(x[1] for x in fwd)
	min_fdG = min(x[2] for x in fwd if x[2] != float('-inf'))
	max_rlength = max(x[1] for x in rev)
	min_rdG = min(x[2] for x in rev if x[2] != float('-inf'))

	pad = 21+max_rlength
	r_fmt = '{: >10s}=ss \"{:->'+str(max_rlength)+'s}\" {:2d}bp'
	m_fmt = ' :{}-{:'+str(math.ceil(math.log(len(fwd),10)))+'d}-{}: '
	f_fmt = '{:2d}bp \"{:->'+str(max_flength)+'s}\" ss={: <10s}'
	for i, ((f_base, f_length, f_dG, f_Tm), 
					(r_base, r_length, r_dG, r_Tm)) in enumerate(zip(fwd,reversed(rev))):
		j = len(sequence)-i
		mval = m_fmt.format(r_base,
												i,
												f_base) 
		if f_length > 0:
			fval = f_fmt.format(f_length, 
													sequence[i-f_length:i], 
													'#'*int(round(10.*f_dG/min_fdG)))
		else:
			fval = '' 
		if r_length > 0:
			rval = r_fmt.format('#'*int(round(10.*r_dG/min_rdG)),
													rsequence[j-r_length:j], 
													r_length)
		else:
			rval = ' '*pad

		print(rval + mval + fval)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Discover potential primer'
																							+'sites in a DNA sequence')

	parser.add_argument('sequence', 
											type=lambda x: is_valid_sequence(parser, x),
											help='The sequence to search')

	parser.add_argument('-Tm',
										 type=float,
										 default=62.0,
										 help='Target primer melting temperature, '+
													'used to calculate length')

	parser.add_argument('-Na',
										 type=float,
										 default=50.0,
										 help='Na concentration for melting temperature')
	parser.add_argument('-K',
										 type=float,
										 default=0.0,
										 help='K concentration for melting temperature')
	parser.add_argument('-Tris',
										 type=float,
										 default=0.0,
										 help='Tris concentration for melting temperature')
	parser.add_argument('-Mg',
										 type=float,
										 default=2.0,
										 help='Mg++ concentration for melting temperature')
	parser.add_argument('-dNTPs',
										 type=float,
										 default=0.2,
										 help='Na concentration for melting temperature')
	parser.add_argument('-oligo',
										 type=float,
										 default=500.,
										 help='Oligo concentration for melting temperature')

	args = parser.parse_args()

	kwargs = vars(args)
	sequence = kwargs.pop('sequence')
	print_output(sequence,
							 analyse(sequence, **kwargs),
							 analyse(reverse_complement(sequence), **kwargs))
