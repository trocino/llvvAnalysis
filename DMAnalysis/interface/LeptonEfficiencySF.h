#ifndef LeptonEfficiencySF_h
#define LeptonEfficiencySF_h

class LeptonEfficiencySF {
public:
  LeptonEfficiencySF(int era=2015) : era_(era) { }
  ~LeptonEfficiencySF() {}

  // pair<SF, pair<stat(SF), syst(SF)> >

  // ================================================================================================

  std::pair<float, std::pair<float,float> > getLeptonEfficiency(float pt, float eta, int id) {

    std::pair<float, std::pair<float,float> > eff(1.0, std::pair<float,float>(0.01, 0.01));

    if(pt<20) return eff; 

    if(id==11) return getLeptonIdEfficiency(pt, eta, id); 

    if(era_==2015) { 
      std::pair<float, std::pair<float,float> > effId  = getLeptonIdEfficiency(pt, eta, id); 
      std::pair<float, std::pair<float,float> > effIso = getLeptonIsoEfficiency(pt, eta, id); 

      eff.first = effId.first * effIso.first; 

      // N.B.: assuming that the errors on ID and ISO scale factors are fully uncorrelated --> certainly not true!! 
      //       more conservative approach would be to treat them as fully correlated (i.e. sum them linearly) 
      eff.second.first  = eff.first * sqrt( (effId.second.first *effId.second.first )/(effId.first*effId.first) + (effIso.second.first *effIso.second.first )/(effIso.first*effIso.first) ); 
      eff.second.second = eff.first * sqrt( (effId.second.second*effId.second.second)/(effId.first*effId.first) + (effIso.second.second*effIso.second.second)/(effIso.first*effIso.first) ); 
    } 

    return eff; 
  }



  // ================================================================================================

  std::pair<float, std::pair<float,float> > getLeptonIdEfficiency(float pt, float eta, int id) {

    std::pair<float, std::pair<float,float> > eff(1.0, std::pair<float,float>(0.01, 0.01));

    if(pt<20) return eff; 

    if(era_==2015) { 

      // Medium ele - ID efficiency scale factors 
      if(id==11) { 
        if(fabs(eta)<0.8) { 
          if     (pt<20) { eff.first = 1.007259; eff.second.first = 0.0656124; eff.second.second = 0.0100000; } 
          else if(pt<30) { eff.first = 0.971182; eff.second.first = 0.0279776; eff.second.second = 0.0100000; } 
          else if(pt<40) { eff.first = 0.984906; eff.second.first = 0.0291216; eff.second.second = 0.0100000; } 
          else if(pt<50) { eff.first = 0.985899; eff.second.first = 0.0175478; eff.second.second = 0.0100000; } 
          else           { eff.first = 0.988598; eff.second.first = 0.0060336; eff.second.second = 0.0100000; } 
        } 
        else if(fabs(eta)<1.444) { 
          if     (pt<20) { eff.first = 1.090196; eff.second.first = 0.0391175; eff.second.second = 0.0100000; } 
          else if(pt<30) { eff.first = 0.983359; eff.second.first = 0.0256741; eff.second.second = 0.0100000; } 
          else if(pt<40) { eff.first = 0.987179; eff.second.first = 0.0171527; eff.second.second = 0.0100000; } 
          else if(pt<50) { eff.first = 0.986920; eff.second.first = 0.0277844; eff.second.second = 0.0100000; } 
          else           { eff.first = 0.986159; eff.second.first = 0.0054099; eff.second.second = 0.0100000; } 
        } 
        else if(fabs(eta)<1.566) { 
          if     (pt<20) { eff.first = 1.086420; eff.second.first = 0.0800687; eff.second.second = 0.0100000; } 
          else if(pt<30) { eff.first = 0.963054; eff.second.first = 0.2618842; eff.second.second = 0.0100000; } 
          else if(pt<40) { eff.first = 0.949123; eff.second.first = 0.0511487; eff.second.second = 0.0100000; } 
          else if(pt<50) { eff.first = 0.981612; eff.second.first = 0.0362270; eff.second.second = 0.0100000; } 
          else           { eff.first = 0.997257; eff.second.first = 0.1134655; eff.second.second = 0.0100000; } 
        } 
        else if(fabs(eta)<2) { 
          if     (pt<20) { eff.first = 0.984444; eff.second.first = 0.0447765; eff.second.second = 0.0100000; } 
          else if(pt<30) { eff.first = 0.936809; eff.second.first = 0.0186922; eff.second.second = 0.0100000; } 
          else if(pt<40) { eff.first = 0.975066; eff.second.first = 0.0274652; eff.second.second = 0.0100000; } 
          else if(pt<50) { eff.first = 0.992806; eff.second.first = 0.0404665; eff.second.second = 0.0100000; } 
          else           { eff.first = 1.005787; eff.second.first = 0.0125727; eff.second.second = 0.0100000; } 
        } 
        else { 
          if     (pt<20) { eff.first = 1.035573; eff.second.first = 0.0303603; eff.second.second = 0.0100000; } 
          else if(pt<30) { eff.first = 0.986446; eff.second.first = 0.0263447; eff.second.second = 0.0100000; } 
          else if(pt<40) { eff.first = 0.963351; eff.second.first = 0.0445603; eff.second.second = 0.0100000; } 
          else if(pt<50) { eff.first = 1.006112; eff.second.first = 0.0103732; eff.second.second = 0.0100000; } 
          else           { eff.first = 1.009490; eff.second.first = 0.0137829; eff.second.second = 0.0100000; } 
        } 
      } // end id==11 

      // Medium mu - ID efficiency scale factors 
      if(id==13) { 
        if(fabs(eta)<0.9) { 
          if     (pt<25) { eff.first = 0.979837; eff.second.first = 0.0027643; eff.second.second = 0.0100000; } 
          else if(pt<30) { eff.first = 0.983967; eff.second.first = 0.0013001; eff.second.second = 0.0100000; } 
          else if(pt<40) { eff.first = 0.989050; eff.second.first = 0.0004463; eff.second.second = 0.0100000; } 
          else if(pt<50) { eff.first = 0.990462; eff.second.first = 0.0003327; eff.second.second = 0.0100000; } 
          else if(pt<60) { eff.first = 0.988163; eff.second.first = 0.0008921; eff.second.second = 0.0100000; } 
          else           { eff.first = 0.989260; eff.second.first = 0.0018805; eff.second.second = 0.0100000; } 
        } 
        else if(fabs(eta)<1.2) { 
          if     (pt<25) { eff.first = 0.991098; eff.second.first = 0.0041236; eff.second.second = 0.0100000; } 
          else if(pt<30) { eff.first = 0.984282; eff.second.first = 0.0023295; eff.second.second = 0.0100000; } 
          else if(pt<40) { eff.first = 0.992365; eff.second.first = 0.0007710; eff.second.second = 0.0100000; } 
          else if(pt<50) { eff.first = 0.992595; eff.second.first = 0.0005348; eff.second.second = 0.0100000; } 
          else if(pt<60) { eff.first = 0.990194; eff.second.first = 0.0014071; eff.second.second = 0.0100000; } 
          else           { eff.first = 0.993889; eff.second.first = 0.0030669; eff.second.second = 0.0100000; } 
        } 
        else if(fabs(eta)<2.1) { 
          if     (pt<25) { eff.first = 0.993387; eff.second.first = 0.0020866; eff.second.second = 0.0100000; } 
          else if(pt<30) { eff.first = 0.992218; eff.second.first = 0.0011782; eff.second.second = 0.0100000; } 
          else if(pt<40) { eff.first = 0.993519; eff.second.first = 0.0004524; eff.second.second = 0.0100000; } 
          else if(pt<50) { eff.first = 0.994435; eff.second.first = 0.0003004; eff.second.second = 0.0100000; } 
          else if(pt<60) { eff.first = 0.991820; eff.second.first = 0.0008790; eff.second.second = 0.0100000; } 
          else           { eff.first = 0.992833; eff.second.first = 0.0022404; eff.second.second = 0.0100000; } 
        } 
        else { 
          if     (pt<25) { eff.first = 0.962504; eff.second.first = 0.0041866; eff.second.second = 0.0100000; } 
          else if(pt<30) { eff.first = 0.963459; eff.second.first = 0.0026978; eff.second.second = 0.0100000; } 
          else if(pt<40) { eff.first = 0.967285; eff.second.first = 0.0013179; eff.second.second = 0.0100000; } 
          else if(pt<50) { eff.first = 0.965518; eff.second.first = 0.0012022; eff.second.second = 0.0100000; } 
          else if(pt<60) { eff.first = 0.960329; eff.second.first = 0.0030173; eff.second.second = 0.0100000; } 
          else           { eff.first = 0.965899; eff.second.first = 0.0078324; eff.second.second = 0.0100000; } 
        } 
      } // end id==13 
    } // end era==2015 

    return eff;
  } 


  // ================================================================================================

  std::pair<float, std::pair<float,float> > getLeptonIsoEfficiency(float pt, float eta, int id) {

    std::pair<float, std::pair<float,float> > eff(1.0, std::pair<float,float>(0.01, 0.01));

    if(pt<20) return eff; 

    if(era_==2015) { 

      // Medium ele - ISO efficiency scale factors 
      if(id==11) { eff.first=1.000; eff.second.first=0.000; eff.second.second=0.000; } // useless, everything done in getLeptonIdEfficiency 


      // Medium mu - ISO efficiency scale factors 
      if(id==13) { 
        if(fabs(eta)<0.9) { 
          if     (pt<25) { eff.first = 1.000248; eff.second.first = 0.0037157; eff.second.second = 0.0050000; } 
          else if(pt<30) { eff.first = 0.997164; eff.second.first = 0.0021178; eff.second.second = 0.0050000; } 
          else if(pt<40) { eff.first = 1.000203; eff.second.first = 0.0007556; eff.second.second = 0.0050000; } 
          else if(pt<50) { eff.first = 0.998514; eff.second.first = 0.0002833; eff.second.second = 0.0050000; } 
          else if(pt<60) { eff.first = 0.998346; eff.second.first = 0.0007213; eff.second.second = 0.0050000; } 
          else           { eff.first = 0.999508; eff.second.first = 0.0008823; eff.second.second = 0.0050000; } 
        } 
        else if(fabs(eta)<1.2) { 
          if     (pt<25) { eff.first = 1.007671; eff.second.first = 0.0057249; eff.second.second = 0.0050000; } 
          else if(pt<30) { eff.first = 0.995432; eff.second.first = 0.0037510; eff.second.second = 0.0050000; } 
          else if(pt<40) { eff.first = 0.999852; eff.second.first = 0.0013890; eff.second.second = 0.0050000; } 
          else if(pt<50) { eff.first = 0.999893; eff.second.first = 0.0004609; eff.second.second = 0.0050000; } 
          else if(pt<60) { eff.first = 1.000990; eff.second.first = 0.0012257; eff.second.second = 0.0050000; } 
          else           { eff.first = 1.000073; eff.second.first = 0.0014886; eff.second.second = 0.0050000; } 
        } 
        else if(fabs(eta)<2.1) { 
          if     (pt<25) { eff.first = 0.993886; eff.second.first = 0.0028343; eff.second.second = 0.0050000; } 
          else if(pt<30) { eff.first = 1.001056; eff.second.first = 0.0018839; eff.second.second = 0.0050000; } 
          else if(pt<40) { eff.first = 1.001161; eff.second.first = 0.0007861; eff.second.second = 0.0050000; } 
          else if(pt<50) { eff.first = 0.998999; eff.second.first = 0.0003977; eff.second.second = 0.0050000; } 
          else if(pt<60) { eff.first = 0.999101; eff.second.first = 0.0006777; eff.second.second = 0.0050000; } 
          else           { eff.first = 1.001630; eff.second.first = 0.0009203; eff.second.second = 0.0050000; } 
        } 
        else { 
          if     (pt<25) { eff.first = 0.994855; eff.second.first = 0.0051735; eff.second.second = 0.0050000; } 
          else if(pt<30) { eff.first = 0.995934; eff.second.first = 0.0034131; eff.second.second = 0.0050000; } 
          else if(pt<40) { eff.first = 0.996796; eff.second.first = 0.0014729; eff.second.second = 0.0050000; } 
          else if(pt<50) { eff.first = 0.999471; eff.second.first = 0.0006253; eff.second.second = 0.0050000; } 
          else if(pt<60) { eff.first = 1.000297; eff.second.first = 0.0014962; eff.second.second = 0.0050000; } 
          else           { eff.first = 0.999578; eff.second.first = 0.0021686; eff.second.second = 0.0050000; } 
        } 
      } // end id==13 
    } // end era==2015 

    return eff;
  } 


private: 
  int era_; 
}; 


#endif
