{
 gROOT->Reset();

 gSystem->Load("libEXOUtilities");

 TChain t("tree");
 t.Add("/nfs/slac/g/exo_data2/exo_data/data/WIPP/root/5851/run00005851-002.root");
 EXOEventData* ed;
 t.SetBranchAddress("EventBranch",&ed);
 t.GetEntry(0);
 EXOWaveformData* wvs = ed->GetWaveformData();
 wvs->Decompress();
 EXOWaveform* wv = wvs->GetWaveformWithChannel(0);
 TH1D* h = wv->GimmeHist();
 h->SetNameTitle("h","Original");
 h->SetLineColor(4);

 
 TChain t2("tree");
 t2.Add("/nfs/slac/g/exo-userdata/users/shaolei/exosim_data/temp_full_5851_real_nocut.root");
 EXOEventData* ed2;
 t2.SetBranchAddress("EventBranch",&ed2);
 t2.GetEntry(0);
 EXOWaveformData* wvs2 = ed2->GetWaveformData();
 wvs2->Decompress();
 EXOWaveform* wv2 = wvs2->GetWaveformWithChannel(0);
 TH1D* h2 = wv2->GimmeHist();
 h2->SetNameTitle("h2","Copy");
 h2->SetLineColor(2);
 
 TCanvas c;
 h->Draw();
 h2->Draw("same");

 c.SaveAs("./test0.png");




}
