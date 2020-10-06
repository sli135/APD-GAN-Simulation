import os, re, sys, glob, math, random
import numpy as np
import ROOT

def main(runNumber, maxEvents,energyThreshold,applyFiducialCut, calibP0, calibP1, outdir):
# LB run limits and file list
    rotationAngle = 0.5312
    time_of_events = 1357754448
    nchannels = 226
    #Create ROOT file and make the Event Data Branch
    standard_wf_length = 2048
    standard_period = 1*ROOT.CLHEP.microsecond
    NoiseFile = os.path.join(outdir, "LB_noise_library_run%i.root" % rn)

    golden_raw_LB_files = []
    golden_masked_LB_files = []

    waveformFiles = '/nfs/slac/g/exo_data[wildcard]/exo_data/data/WIPP/root/[runNumber]/'
    preprocessedFiles = '/nfs/slac/g/exo_data4/users/Energy/data/WIPP/preprocessed/2017_Phase1_v2/ForCalibration/run_[runNumber]_tree.root'
    ROOT.gSystem.Load('/nfs/slac/g/exo/shaolei/EXOEnergy/lib/libEXOEnergy.so')
    infoRuns = ROOT.EXOSourceRunsPolishedInfo('/nfs/slac/g/exo/shaolei/EXOEnergy/data/SourceRunsInfo_Phase1_20170411.txt')
    print('Opening preprocessed file %s'%(preprocessedFiles.replace('[runNumber]','%i'%runNumber)))
    fProcessed = ROOT.TFile(preprocessedFiles.replace('[runNumber]','%i'%runNumber),'READ')
    tProcessed = fProcessed.Get('dataTree')
    
    print('Opening waveform files in directory %s'%(waveformFiles.replace('[runNumber]','%i'%runNumber)))
    tWF = ROOT.TChain('tree','tree')

    for i in range(2,8):
        wfFile = waveformFiles.replace('[wildcard]','%i'%i)+'*.root'
        tWF.Add(wfFile.replace('[runNumber]','%i'%runNumber))
    
    ed_WF = ROOT.EXOEventData()
    tWF.SetBranchAddress('EventBranch',ed_WF)
    
    tWF.BuildIndex("EventBranch.fRunNumber","EventBranch.fEventNumber")

    #cut = "nsc == 1 && multiplicity == 1 && isSolicitedTrigger == 0 && t_scint > -250 && t_scint < 100"
    cut = "nsc == 1 && multiplicity == 1"

    if applyFiducialCut:
        cut += " && TMath::Abs(cluster_z) > 10.0 && TMath::Abs(cluster_z) < 182.0 && TMath::Sqrt(cluster_x*cluster_x+cluster_y*cluster_y) < 162.0"
        #cut += " && (cluster_x > 0 || (cluster_x < 0 && TMath::Abs(cluster_y) > TMath::Abs(cluster_x))) && TMath::Abs(cluster_z) > 10.0"

    if energyThreshold > 0.0:
        cut += " && (e_charge*TMath::Cos(%f)+e_scint*TMath::Sin(%f))*0.734 + 11.8 > %f"%(rotationAngle,rotationAngle,energyThreshold)

    # Diagonal cut
    cut += " && e_scint < e_charge*1.1 + 617.4"

    n = tProcessed.Draw("eventNum:e_scint:t_scint",cut,"goff")
    eventNumbers = np.copy(np.frombuffer(tProcessed.GetV1(), count=n))
    e_scint = np.copy(np.frombuffer(tProcessed.GetV2(), count=n))
    t_scint = np.copy(np.frombuffer(tProcessed.GetV3(), count=n))

    n = tProcessed.Draw("cluster_x:cluster_y:cluster_z:(%f*e_charge+%f)"%(calibP1,calibP0),cut,"goff")
    cluster_x = np.copy(np.frombuffer(tProcessed.GetV1(), count=n))
    cluster_y = np.copy(np.frombuffer(tProcessed.GetV2(), count=n))
    cluster_z = np.copy(np.frombuffer(tProcessed.GetV3(), count=n))
    e_charge = np.copy(np.frombuffer(tProcessed.GetV4(), count=n))
    
    # Build index event number
    tWF.BuildIndex("fRunNumber", "fEventNumber")
    # The following is a list of tuples of (run, event)
    v1 = tProcessed.GetV1()
    v2 = tProcessed.GetV2()
    run_solicited_events = [(v1[i], v2[i]) for i in range(tProcessed.GetSelectedRows())]
    
    # Make a Double WF Array for the noise glitch templates
    #events_array = MakeGlitchTempArray(tWF, run_solicited_events)
    # Create new root file
    my_save_path = "/nfs/slac/g/exo-userdata/users/shaolei/exosim_data/"
    open_file = ROOT.TFile(my_save_path+"temp_full_%i_real.root" % runNumber, "recreate")
    save_tree = ROOT.TTree("tree", "tree")
    ed = ROOT.EXOEventData()
    save_tree.Branch("EventBranch", ed)

    print(save_tree.GetEntries())
    print(tWF.GetEntries())

    count = 0
    k = 0
    start_point = k * 10000
    pickle_dir = '/nfs/slac/g/exo-userdata/users/shaolei/data_pickles/'
    xdata,ydata = [],[]
    nevents = 10000
    count, end_point = 0, start_point + nevents
    with open(pickle_dir+'GAN_generated_waveform_5851.csv','rb') as f:
    #with open('/nfs/slac/g/exo-userdata/users/shaolei/data/WF/EnergyThreshold_800_calib_margin2_v2/APDWFSignalsEnergyThreshold_5851.csv','rb') as f:
      for line in f:
              if count > end_point:
                break
              if count > start_point:
                x = np.fromstring(line, dtype=float, sep=',')
                xdata.append(x[0:74*350])
                ydata.append(x[74*350:])
              count += 1
    ydata = np.array(ydata)
    print(ydata.shape)
    print(n)
    #print(ydata[:10])
    #print('events_array',events_array)
    print('e_charge',e_charge[:10])
    print('e_scint',e_scint[:10])
    print('run_solicited_events',run_solicited_events[:10])
    xdata = np.array(xdata) + 1600
    if nevents > xdata.shape[0]:
        nevents = xdata.shape[0]
    xdata = xdata.reshape(nevents,74,350)
    for i in range(n):
        tWF.GetEntryWithIndex(runNumber,int(eventNumbers[i]))
        raw_ED = tWF.EventBranch
        wf_data = raw_ED.GetWaveformData()
        wf_data.Decompress()
        if int(wf_data.GetNumWaveforms()) != nchannels: 
            #All 226 channel need to exist
            continue

        if i % 100 == 0: print('%i events processed.' % i)
        ed.fRunNumber = runNumber
        ed.fEventNumber = i
        ed.fEventHeader.fTriggerSeconds = time_of_events
        for chan in range(226):
           
            #Create EXO Waveform from numpy array
            
            if chan in []:#range(152,226):
              data = xdata[i][chan - 152]
              arr = data
              mu, sigma = np.mean(arr[:50]), np.std(arr[:50])
              noise_1 = np.random.normal(mu, sigma,size=750)
              noise_2 = np.random.normal(mu, sigma,size=948)
              arr = np.concatenate([noise_1,data])
              arr = np.concatenate([arr,noise_2])
              #arr = arr.astype(np.int64)
              #print(chan,arr.dtype)
            else:
              waveformData = wf_data.GetWaveform(chan)
              if not waveformData:
                print('waveform does not exist.')
              arr = np.array(waveformData)
              arr = arr.astype(np.float)
              #print(chan,arr.dtype)
            #print(chan,arr,arr.shape)
            dw_wf = ROOT.EXODoubleWaveform(arr, len(arr))
            dw_wf.SetSamplingPeriod(standard_period)
            exo_wf = ed.GetWaveformData().GetNewWaveform()
            exo_wf.ConvertFrom(dw_wf)
            exo_wf.fChannel = chan
            dw_wf.IsA().Destructor( dw_wf )
        ed.GetWaveformData().fNumSamples = standard_wf_length
        ed.GetWaveformData().Compress()
        save_tree.Fill() #Fill tree with ED
        #ed.Clear() #Clear ED for next event
    if not save_tree:
        print('Did not create tree successfully.')
    else:
        print('tree Entries:',save_tree.GetEntries())
    save_tree.Write() #Write to the file
    open_file.Close() 
    open_file.IsA().Destructor( open_file )
    #Loop over all channels (Charge + Light)

# Run script standalone
if __name__ == "__main__":

    if ROOT.gSystem.Load("libEXOROOT") < 0: sys.exit('Failed to load EXOROOT.')
    if ROOT.gSystem.Load("libEXOUtilities") < 0: sys.exit('Failed to load EXOUtilities')
  
    if ROOT.gSystem.Load('../EXOEnergy/lib/libEXOEnergy.so') < 0: sys.exit('Failed to load libEXOEnergy.so')
  
    if len(sys.argv) < 4:
        print "arguments: [run number] [events to use from run] [energy threshold] [out dir]"
        sys.exit(1)
    ###################### Input parameters ########################
    rn = int(sys.argv[1])
    nevents = int(sys.argv[2])
    energyThreshold = float(sys.argv[3])
    outdir = sys.argv[4]
    ################################################################
    infoRuns = ROOT.EXOSourceRunsPolishedInfo('../EXOEnergy/data/SourceRunsInfo_Phase1_20170411.txt')
    
    allRuns = infoRuns.GetListOf('RunNumber')
    #runNumbers = [int(run) for run in allRuns if int(run) < 3000]

    runs = set([])
    for i in range(allRuns.size()):
        runs.add(allRuns.at(i))
    #runNumbers = [int(run) for run in runs if int(run) < 5522 and int(run) > 5422]
    calib = ROOT.TTree()
    calib.ReadFile('phase1_calibPars_ionization.txt')
    applyCalibration = True
    
    info = infoRuns.Clone()
    info.CutExact('RunNumber','%s'%rn)
    week = info.GetListOf('WeekIndex')[0]
    calib.Draw('p0:p1','week == %s && multiplicity == 1'%week,'goff')
    
    calibP0 = 0.0
    calibP1 = 1.0
    if applyCalibration:
        calibP0 = calib.GetV1()[0]
        calibP1 = calib.GetV2()[0]
    print('calibP0',calibP0,'calibP1',calibP1)
    applyFiducialCut = True
    print "Process %i events for  run %i " % (nevents, rn)
    main(rn, nevents, energyThreshold, applyFiducialCut, calibP0, calibP1, outdir)

