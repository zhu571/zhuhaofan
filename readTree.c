//将下列代码保存到readTree.cc
//在ROOT环境下 .x readTree.cc  

void readTree()
{
// 1.打开文件，得到TTree指针
  TFile *ipf=new TFile("tree.root");//打开ROOT文件
  if (ipf->IsZombie()) {
   cout << "Error opening file" << endl;
   exit(-1);
  }
  ipf->cd();
  TTree *tree=(TTree*)ipf->Get("tree");//得到名字为“tree”的TTree指针

//2. 声明tree的Branch变量
  Double_t x;
  Double_t e;
  int pid;
  Double_t tof, ctof;
  Double_t tu, td;
  Double_t qu, qd;

//3. 将变量指向对应Branch的地址
  tree->SetBranchAddress("ctof",&ctof);//将ROOT文件内tree内名为"ctof"的branch的数据的指针指向ctof的变量。
  tree->SetBranchAddress("tof",&tof);  
  tree->SetBranchAddress("pid",&pid);
  tree->SetBranchAddress("tu",&tu);   
  tree->SetBranchAddress("td",&td);
  tree->SetBranchAddress("qu",&qu);   
  tree->SetBranchAddress("qd",&qd);
  tree->SetBranchAddress("x",&x);

  //Histogram
  TH1D *hTOF= new TH1D("hTOF","Time of flight", 1000,0,100);
  TH1D *tdiff=new TH1D("tdiff","td-tu",140,-20,50);  
  TH1D *qdiff=new TH1D("qdiff","log(qu/qd)",200,-1,1);
  TH1D *dtd=new TH1D("dtd","dt/dx",141,-20.25,50.25);
  TH1D *dqd=new TH1D("dqd","dt/dx",100,-1,1);
  TH1D *htx=new TH1D("htx","htx",500,-120,120);
  TH2D *hdx=new TH2D("hdx","htx-hx:hx",100,-20,20,500,-120,120);
  TH1D *qhtx=new TH1D("qhtx","qhtx",500,-120,120);
  TH2D *qhdx=new TH2D("qhdx","qhtx-hx:hx",400,-80,80,500,-120,120);
  TH2D *hgtofx=new TH2D("hgtofx","hgtofx",100,-120,120,100,39,48);
  TH1D *hgctof=new TH1D("hgctof","hgctof",100,39,48);
  TH2D *hgtofcx=new TH2D("hgtofcx","corrected TOF",100,-120,120,100,15,19);
  TH1D *htofc=new TH1D("htofc","htofc",200,0,100);//calculated tof - normalized to 500cm
  TH1D *htof=new TH1D("htof","htof",200,0,100);//real tof - normalized to 500cm
  //qx修正
  TH2D *qhgtofx=new TH2D("qhgtofx","qhgtofx",100,-120,120,100,39,48);
  TH1D *qhgctof=new TH1D("qhgctof","qhgctof",100,39,48);
  TH2D *qhgtofcx=new TH2D("qhgtofcx","corrected TOF",100,-120,120,100,15,19);
  TH1D *qhtofc=new TH1D("qhtofc","qhtofc",200,0,100);//calculated tof - normalized to 500cm
  TH1D *qhtof=new TH1D("qhtof","qhtof",200,0,100);//real tof - normalized to 500cm

  //将新数据写入新的ROOT文件 -对应的代码用 ////标出
  //// //calibration parameters
  //// Double_t a,b;
  //// ... ... ...
  //// //new tree parameters
  Double_t tx,qx,ntof,ce;
  //// redefinition... ... ...    
  TFile *opf=new TFile("tree2.root","recreate");
  opf->cd();
  TTree *opt=new TTree("tree","tree");
  opt->Branch("tx",&tx,"tx/D");
  opt->Branch("qx",&qx,"qx/D");
  opt->Branch("ntof",&ntof,"ntof/D");
  opt->Branch("ce",&ce,"ce/D");
  //// ... ... ...


  //4. 逐事件读取tree的branch数据
  Long64_t nentries=tree->GetEntries();//得到tree的事件总数
  for(Long64_t jentry=0; jentry<nentries; jentry++) {//对事件进行遍历
    tree->GetEntry(jentry);//将第jentry个事件数据填入对应变量(步骤3.中指向的变量)，每次变量值会变成当前事件对应的数据。
    hTOF->Fill(ctof);
    tdiff->Fill(td-tu);
    qdiff->Fill(log(qu/qd));
    //// // calculate new parameters
    tx=tu-td;
    qx=log(qu/qd);
    //// ... ... ...

    //ntof
    Double_t x=3.742*(td-tu-14.93);
    Double_t d=TMath::Sqrt(502.5*502.5+x*x);

    ntof=(ctof-26.18)/d*100.;
    //中子能量
    if(ctof>45) { //选择中子
	ce=72.29824*72.29824/(ntof*ntof);

    }




    opt->Fill();//fill new parameter to TTree* opt

    if(jentry%100000==0) cout<<"process "<<jentry<<" of "<<nentries<<endl;
  }
  //tx导数
    for(int i=1;i<=tdiff->GetNbinsX();i++) {
         Double_t df=abs(tdiff->GetBinContent(i+1)-tdiff->GetBinContent(i));
         dtd->Fill(tdiff->GetBinLowEdge(i+1),df);
  }
    //tx修正
    for(Long64_t jentry=0; jentry<nentries; jentry++) {//对事件进行遍历
    tree->GetEntry(jentry);
    Double_t tx=3.750*(td-tu-14.90);
    htx->Fill(tx);
    hdx->Fill(tx-x,x);//difference
  }
    //qx导数
    for(int i=1;i<=qdiff->GetNbinsX();i++) {
          Double_t qdf=abs(qdiff->GetBinContent(i+1)-qdiff->GetBinContent(i));
           dqd->Fill(qdiff->GetBinLowEdge(i+1),qdf);
  }

    //qx修正
    for(Long64_t jentry=0; jentry<nentries; jentry++) {//对事件进行遍历
         tree->GetEntry(jentry);
         Double_t qx=189.7*log(qu/qd);
         qhtx->Fill(qx);
         qhdx->Fill(qx-x,x);//difference
   }
    //关联一维图
  TH1D *hdx1=hdx->ProjectionX("projx of hdx");
  TH1D *qhdx1=qhdx->ProjectionX("projx of qhdx");
   
   //TOF修正
   for(Long64_t jentry=0; jentry<nentries; jentry++) {//对事件进行遍历
    tree->GetEntry(jentry);
    Double_t tx=3.742*(td-tu-14.93);
    if(ctof>42&& ctof<44.5) { //选择gamma
        Double_t d=TMath::Sqrt(502.5*502.5+tx*tx);
        Double_t ctofa=(ctof)/d*500.;//normalized to 500cm
        hgtofx->Fill(tx,ctofa);
        if(abs(tx)<5) hgctof->Fill(ctof);//gamma hits the center of the det.
    }
  }
  //TOF绝对刻度
  for(Long64_t jentry=0; jentry<nentries; jentry++) {
    tree->GetEntry(jentry);
    Double_t tx=3.742*(td-tu-14.93);
    Double_t d=TMath::Sqrt(502.5*502.5+tx*tx);
    Double_t tofc=(ctof-26.18)/d*500.;//normalized to 500cm
    hgtofcx->Fill(tx,tofc);//gamma hits the center of the det.
    htofc->Fill(tofc);
    htof->Fill(tof*500./d);
  } 

  //TOF qx修正
  for(Long64_t jentry=0; jentry<nentries; jentry++) {//对事件进行遍历
    tree->GetEntry(jentry);
    Double_t qx=189.7*(log(qu/qd));
    if(ctof>42&& ctof<44.5) { //选择gamma
        Double_t qd=TMath::Sqrt(502.5*502.5+qx*qx);
        Double_t qctofa=(ctof)/qd*500.;//normalized to 500cm
        qhgtofx->Fill(qx,qctofa);
        if(abs(qx)<5) qhgctof->Fill(ctof);//gamma hits the center of the det.
    }
  }
  //TOF qx绝对刻度
  for(Long64_t jentry=0; jentry<nentries; jentry++) {
    tree->GetEntry(jentry);
    Double_t qx=189.7*(log(qu/qd));
    Double_t qd=TMath::Sqrt(502.5*502.5+qx*qx);
    Double_t qtofc=(ctof-26.18)/qd*500.;//normalized to 500cm
    qhgtofcx->Fill(qx,qtofc);//gamma hits the center of the det.
    qhtofc->Fill(qtofc);
    qhtof->Fill(tof*500./qd);
  }


  //结束部分
  hTOF->Write();
  hdx1->Write();
  hdx1->Fit("gaus");
  qhdx1->Write();
  qhdx1->Fit("gaus");
  htx->Write();
  hdx->Write();
  qhtx->Write();
  qhdx->Write();
  dtd->Sumw2(0);//不显示传递误差
  dtd->Write();
  dqd->Sumw2(0);//不显示传递误差
  dqd->Write();
  dqd->Fit("gaus","","",-0.7,-0.35);
  dqd->Fit("gaus","","",0.3,0.7);
  tdiff->Write();
  qdiff->Write();
  hgtofx->Write();
  hgctof->Write();
  hgctof->Fit("gaus");
  hgtofcx->Write();
  htofc->Write();
  htof->Write();
  qhgtofx->Write();
  qhgctof->Write();
  qhgctof->Fit("gaus");
  qhgtofcx->Write();
  qhtofc->Write();
  qhtof->Write();
  ipf->Close();
  opt->Write();
  opf->Close();
}
