#!/usr/bin/python
# vim: set encoding=utf-8

# usage: python $0 samples.csv | cat header.csv - > sample_sheet.csv

import pandas as pd
import sys
import string

df = pd.read_csv(sys.argv[1], header=0)

samples = df.T.iloc[1:,:]
samples.columns = ['run', 'plate', 'array', 'column', 'plate_column'] + map(str, range(1, 7))
samples = pd.melt(samples, id_vars=['run', 'plate', 'array', 'column', 'plate_column'], var_name='row', value_name='sample_id')

samples.run = samples.run.astype('int')
samples.row = samples.row.astype('int')
samples.column = samples.column.astype('int')
samples.plate_column = samples.plate_column.astype('int')
samples.sample_id = samples.sample_id.map(lambda s: s.strip().upper().replace(' ', ''))

samples['unique_id'] = ''
for sample, count in samples.sample_id.value_counts().iteritems():
	if count > 1:
		samples.loc[samples.sample_id == sample, 'unique_id'] = map(lambda x, y: '{}_{}'.format(x,y), samples[samples.sample_id == sample].sample_id, range(0, count))
	else:
		samples.loc[samples.sample_id == sample, 'unique_id'] = sample

for idx, row in samples.sort(['run', 'array', 'row', 'column']).iterrows():
	d = row.to_dict()
	d['well'] = '{0}{1:02d}'.format(string.letters[d['row']-1].upper(), d['plate_column'])
	print '{unique_id},{sample_id},{run}_{plate},{well},{array},R{row:02d}C{column:02d},,,,,'.format(**d)
