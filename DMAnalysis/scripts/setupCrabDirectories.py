import shutil
import os

inputpath = ''
campaignInfo = 'RunIISpring16MiniAODv1-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v2/MINIAODSIM'
toppath = ''
datasets = []

### Get datasets from input list
with open(inputpath,'r') as f:
	for l in f.readlines():
		if(len(l)):
			datasets.append( tuple(l.split(':')) )

### For each dataset
###		create a separate subfolder
###		put the template config in it
###		fill the template
for ds in datasets:
	thispath = os.path.join( toppath, dtag )
	os.mkdir( thispath )
	shutil.copy( template, os.path.join( thispath, 'crab_' + ds[0] + '.py' ) )
