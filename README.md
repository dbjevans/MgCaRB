MgCaRB: Compute salinity and pH adjusted planktonic Mg/Ca SST

Very quick start guide. Type: help MgCaRB
at the command line for further details.

1. Navigate Matlab path to unzipped MgCaRB_v1 folder

2. At the command line:

data = csvread('testdata_pCO2.csv',1,0);	% load example data
age = data(:,1);
MgCa = data(:,2);
[TOut,pHOut,TOutRel] = MgCaRB_v1(age,MgCa,1,1,[0 10],2350,0,0,35,500,1);	% run function
   

If this function is useful to your research, please cite:
  Gray, W. & Evans, D. [2019] Nonthermal influences on Mg/Ca in planktonic foraminifera: 
    A review of culture studies and application to the last glacial maximum
    Paleoceanography & Paleoclimatology vol.34 doi.org/10.1029/2018PA003517
