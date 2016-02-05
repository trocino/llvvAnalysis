#ifndef LeptonEfficiencySF_h
#define LeptonEfficiencySF_h

class LeptonEfficiencySF {
public:
  LeptonEfficiencySF(int era=2015) : era_(era) { }
  ~LeptonEfficiencySF() {}

  // pair<SF, pair<stat(SF), syst(SF)> >

  std::pair<float, std::pair<float,float> > getLeptonEfficiency(float pt, float eta, int id) {

    std::pair<float, std::pair<float,float> > eff(1.0, std::pair<float,float>(0.01, 0.01));

    if(era_==2015) { 

      if(pt<20) return eff; 

      // Medium ele - ID + ISO efficiency scale factors 
      if(id==11) { eff.first=1.000; eff.second.first=0.000; eff.second.second=0.010; }

      //// When we finally have all efficiencies and errors 
      // if(id==11) {
      // 	if(fabs(eta)<0.8) {
      // 	  if(pt<30)      { eff.first=1.010; eff.second.first=0.000; eff.second.second=0.010; }
      // 	  else if(pt<40) { eff.first=1.006; eff.second.first=0.000; eff.second.second=0.010; }
      // 	  else if(pt<50) { eff.first=1.009; eff.second.first=0.000; eff.second.second=0.010; }
      // 	  else           { eff.first=1.008; eff.second.first=0.000; eff.second.second=0.010; }
      // 	}
      // 	else if(fabs(eta)<1.442) {
      // 	  if(pt<30)      { eff.first=0.981; eff.second.first=0.000; eff.second.second=0.010; }
      // 	  else if(pt<40) { eff.first=0.987; eff.second.first=0.000; eff.second.second=0.010; }
      // 	  else if(pt<50) { eff.first=0.993; eff.second.first=0.000; eff.second.second=0.010; }
      // 	  else           { eff.first=0.995; eff.second.first=0.000; eff.second.second=0.010; }
      // 	} 
      // 	else if(fabs(eta)<1.556) {
      // 	  if(pt<30)      { eff.first=1.046; eff.second.first=0.000; eff.second.second=0.010; }
      // 	  else if(pt<40) { eff.first=1.011; eff.second.first=0.000; eff.second.second=0.010; }
      // 	  else if(pt<50) { eff.first=0.994; eff.second.first=0.000; eff.second.second=0.010; }
      // 	  else           { eff.first=0.997; eff.second.first=0.000; eff.second.second=0.010; }
      // 	}
      // 	else if(fabs(eta)<2.0) {
      // 	  if(pt<30)      { eff.first=0.992; eff.second.first=0.000; eff.second.second=0.010; }
      // 	  else if(pt<40) { eff.first=0.993; eff.second.first=0.000; eff.second.second=0.010; }
      // 	  else if(pt<50) { eff.first=1.008; eff.second.first=0.000; eff.second.second=0.010; }
      // 	  else           { eff.first=1.009; eff.second.first=0.000; eff.second.second=0.010; }
      // 	}
      // 	else {
      // 	  if(pt<30)      { eff.first=1.045; eff.second.first=0.000; eff.second.second=0.010; }
      // 	  else if(pt<40) { eff.first=1.031; eff.second.first=0.000; eff.second.second=0.010; }
      // 	  else if(pt<50) { eff.first=1.019; eff.second.first=0.000; eff.second.second=0.010; }
      // 	  else           { eff.first=1.014; eff.second.first=0.000; eff.second.second=0.010; }
      // 	}
      // }

      // Medium mu - ID + ISO efficiency scale factors 
      if(fabs(eta)<0.9) { 
	if     (pt<25) { eff.first = 0.979414; eff.second.first = 0.0027415; eff.second.second = 0.0100000; } 
	else if(pt<30) { eff.first = 0.984815; eff.second.first = 0.0012909; eff.second.second = 0.0100000; } 
	else if(pt<40) { eff.first = 0.989551; eff.second.first = 0.0004332; eff.second.second = 0.0100000; } 
	else if(pt<50) { eff.first = 0.991315; eff.second.first = 0.0003197; eff.second.second = 0.0100000; } 
	else if(pt<60) { eff.first = 0.987545; eff.second.first = 0.0008304; eff.second.second = 0.0100000; } 
	else           { eff.first = 0.991321; eff.second.first = 0.0016589; eff.second.second = 0.0100000; } 
      } 
      else if(fabs(eta)<1.2) { 
	if     (pt<25) { eff.first = 0.987974; eff.second.first = 0.0040435; eff.second.second = 0.0100000; } 
	else if(pt<30) { eff.first = 0.985705; eff.second.first = 0.0022352; eff.second.second = 0.0100000; } 
	else if(pt<40) { eff.first = 0.991785; eff.second.first = 0.0007619; eff.second.second = 0.0100000; } 
	else if(pt<50) { eff.first = 0.992362; eff.second.first = 0.0005198; eff.second.second = 0.0100000; } 
	else if(pt<60) { eff.first = 0.991781; eff.second.first = 0.0013460; eff.second.second = 0.0100000; } 
	else           { eff.first = 0.992222; eff.second.first = 0.0028402; eff.second.second = 0.0100000; } 
      } 
      else if(fabs(eta)<2.1) { 
	if     (pt<25) { eff.first = 0.995794; eff.second.first = 0.0021834; eff.second.second = 0.0100000; } 
	else if(pt<30) { eff.first = 0.991305; eff.second.first = 0.0012558; eff.second.second = 0.0100000; } 
	else if(pt<40) { eff.first = 0.993265; eff.second.first = 0.0004708; eff.second.second = 0.0100000; } 
	else if(pt<50) { eff.first = 0.994487; eff.second.first = 0.0003082; eff.second.second = 0.0100000; } 
	else if(pt<60) { eff.first = 0.991223; eff.second.first = 0.0009382; eff.second.second = 0.0100000; } 
	else           { eff.first = 0.995587; eff.second.first = 0.0024502; eff.second.second = 0.0100000; } 
      } 
      else { 
	if     (pt<25) { eff.first = 0.976729; eff.second.first = 0.0047255; eff.second.second = 0.0100000; } 
	else if(pt<30) { eff.first = 0.967756; eff.second.first = 0.0030171; eff.second.second = 0.0100000; } 
	else if(pt<40) { eff.first = 0.967950; eff.second.first = 0.0014332; eff.second.second = 0.0100000; } 
	else if(pt<50) { eff.first = 0.965109; eff.second.first = 0.0013090; eff.second.second = 0.0100000; } 
	else if(pt<60) { eff.first = 0.957879; eff.second.first = 0.0033303; eff.second.second = 0.0100000; } 
	else           { eff.first = 0.974770; eff.second.first = 0.0106513; eff.second.second = 0.0100000; } 
      } 
    } 

    return eff;
  } 


private: 
  int era_; 
}; 


#endif
