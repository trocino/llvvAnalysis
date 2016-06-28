
from WMCore.Configuration import Configuration 

config = Configuration()

config.section_("General")
config.General.requestName = '@NAME'
config.General.workArea = '@WORKAREA'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'ANALYSIS'
config.JobType.psetName = '@CONFIGFILE'

outFile = config.General.requestName.split('_ext')[0]
config.JobType.pyCfgParams = ["outputFile="+outFile]

config.section_("Data")
config.Data.inputDataset = '@INPUTDATASET'
config.Data.lumiMask ='@JSON'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.splitting = '@Method'
config.Data.unitsPerJob = @FILESPERJOB
config.Data.outLFNDirBase = '@OUTPATH'       # e.g. /store/user/aalbert/myFolder 
config.Data.outputDatasetTag = '@OUTTAG'     # e.g. 76X_mcRun2_asymptotic_RunIIFall15DR76_v1

config.section_("Site")
config.Site.storageSite = '@STORAGESITE'

config.section_("User")
#config.User.voGroup = 'dcms'
