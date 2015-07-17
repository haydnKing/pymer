"""Simple python interface to RNAfold"""

import subprocess, re

exp = re.compile(r'\(\s*(-?\d+\.\d+)\)')

def _fold(seq):
	seq = seq.upper().replace('T', 'U')
	s = subprocess.check_output(['RNAfold', '--noPS',],
															input="{}\n@\n".format(seq),
															universal_newlines=True)

	try:
		m = exp.search(s.split('\n')[1])
		f = float(m.groups(1)[0])
	except AttributeError as e:
		raise RuntimeError('Couldn\'t parse RNAfold output: \"{}\"'.format(s))
	return f

def fold(seq):
	if type(seq) == str:
		return _fold(seq)
	else:
		return _fold(str(seq.seq))

