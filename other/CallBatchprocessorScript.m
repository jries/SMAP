SMAP %you need to start SMAP to add directories and set some variables
smapobject=g;

%main batch file:
batchfile='/Volumes/Ries_Data/DataSSD/workingdirectory/GlobalFitPaper/Exampledata/NPC4C/191022_beads/003_beads_10nm_3D_2c_e_676_685_Z-stack_1_batch.mat';
%list of data files:
rawfiles={'/Volumes/Ries_Data/DataSSD/workingdirectory/GlobalFitPaper/Exampledata/NPC4C/191022_beads/001_beads_10nm_3D_2c_e_676_685_Z-stack_1/001_beads_10nm_3D_2c_e_676_685_Z-stack_1_MMStack_Pos0.ome.tif'};

%make batch processor programmatically
bp=WorkflowModules.Batchprocessor;
bp.attachPar(smapobject.P)
bp.makeGui;
bp.guihandles.mainbatchfile.String=batchfile;
bp.guihandles.filelist.String=rawfiles;
bp.guihandles.filelist.Value=1;

%run batch processor:
bp.processb_callback(0,0);