#!/usr/bin/python
# vim: set encoding=utf-8

# usage: python $0 samples.csv | cat header.csv - > sample_sheet.csv

import pandas as pd
import sys
import string

df = pd.read_csv(sys.argv[1], header=None)

samples = df.T.iloc[1:,:]
samples.columns = ['plate', 'array', 'array_column', 'plate_column'] + map(str, range(1, 7))
samples = pd.melt(samples, id_vars=['plate', 'array', 'array_column', 'plate_column'], var_name='row', value_name='sample_id')

# Drop rows with empty sample_id
samples = samples[pd.notnull(samples.sample_id)]

samples.plate = samples.plate.astype('int')
samples.row = samples.row.astype('int')
samples.array_column = samples.array_column.astype('int')
samples.plate_column = samples.plate_column.astype('int')
samples.sample_id = samples.sample_id.map(lambda s: s.strip().upper().replace(' ', ''))

def handle_replicates(replicates):
	rows, cols = replicates.shape
	if rows > 1:
		replicates.unique_id = replicates.sample_id.str.cat(map(str, range(0, rows)), sep='_')
		replicates.replicate = replicates.sample_id.str.cat(('0',)*rows, sep='_')
	return replicates

samples['unique_id'] = samples.sample_id
samples['replicate'] = ''
samples = samples.groupby('sample_id').apply(handle_replicates)
samples.loc[samples.unique_id.str.endswith('_0'),'replicate'] = ''

for idx, row in samples.sort(['plate', 'array', 'row', 'array_column']).iterrows():
	d = row.to_dict()
	d['well'] = '{0}{1:02d}'.format(string.letters[d['row']-1].upper(), d['plate_column'])
	print '{unique_id},{sample_id},{plate},{well},{array},R{row:02d}C{array_column:02d},,,{replicate},,'.format(**d)
