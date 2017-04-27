#include "SHERPA/Tools/Analysis_Interface.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include <algorithm>

/*compile with:
  SHERPA_PREFIX=/home/s0118321/software/Sherpa/rel-2-2-2
  g++ -shared -g -I`$SHERPA_PREFIX/bin/Sherpa-config --incdir`  `$SHERPA_PREFIX/bin/Sherpa-config --ldflags`  -fPIC -o libSELANA.so SELANA.C
 */

namespace SELAN{


  class edr {
  public:
    double E;
    double dr;
    edr(double _e,double _dr) : E(_e), dr(_dr) {}
  };
  class Order_edr {
  public:
    int operator()(const edr a, const edr b);
  };
  int Order_edr::operator()(const edr a, const edr b) {
    if (a.dr<b.dr) return 1;
    return 0;
  }

  class SELANA: public SHERPA::Analysis_Interface {

  protected:

  private:
    bool m_direct, m_check_cuts;
    double m_dr_gamma_lep, m_pt, m_eps, m_n, m_dr;
    std::string m_inpath;
    std::string m_infile;
    std::string m_outpath;


    /*structure of the cut functions:
        return true, if condition for ME-case is fullfilled
        return false, if not. This defines the YFS phase space
    */

    bool LeptonCut(const ATOOLS::Particle & photon, const ATOOLS::Particle_Vector & leptonen){
      /*ME condition: all leptons have to be separated from the hardest photon by m_dr_gamma_lep*/
      for(ATOOLS::Particle_Vector::const_iterator it=leptonen.begin(); it!=leptonen.end(); it++){
           if( photon.Momentum().DR((*it)->Momentum()) < m_dr_gamma_lep) return false;
        }
      return true;
    }

    bool PhotonPtCut(const ATOOLS::Particle &photon){
      /*ME condition: pT(gamma) > m_pt*/
      return (photon.Momentum().PPerp() > m_pt);
    }

    bool Chi(const double & egamma, const double & etot, const double &dr){
      double e_chi = egamma*m_eps*std::pow((1.-std::cos(dr))/(1.-std::cos(m_dr)),m_n);
      if (etot > e_chi) return false;
    return true;
    }

    bool IsoCut(const ATOOLS::Particle &photon, const ATOOLS::Particle_Vector &partonen){
      /*ME condition: etot < e_chi for every radius in the cone */
      const double egamma = photon.Momentum().PPerp();
      std::vector<edr> edrlist;
      for (size_t i=0; i<partonen.size(); i++) {
          ATOOLS::Particle *parton(partonen.at(i));
          double dr =  photon.Momentum().DR(parton->Momentum());
          if (dr<m_dr) edrlist.push_back(edr(parton->Momentum().PPerp(), dr));
      }
      if (!edrlist.empty()){
         std::stable_sort(edrlist.begin(), edrlist.end(), Order_edr());
         double etot=0;
         for (size_t i=0; i< edrlist.size();i++){
             etot+=edrlist[i].E;
             if (!Chi(egamma, etot, edrlist[i].dr)) return false;
           }
        }
    return true;
    }


  public:


    SELANA(const std::string &inpath,
           const std::string &infile,
           const std::string &outpath):SHERPA::Analysis_Interface("SELANA"),
      m_inpath(inpath), m_infile(infile), m_outpath(outpath)
    {
      msg_Debugging()<<"SELANA Selector activ."<<std::endl;
    }

    void ShowSyntax(const int i)
    {
      if (!msg_LevelIsInfo() || i==0) return;
      msg_Out()<<METHOD<<"(): {\n\n"
              <<"custom analysis acting as selector"
             <<"}"<<std::endl;
    }


    bool Init()
    {

      ATOOLS::Data_Reader reader(" ",";","//","=");
      reader.AddWordSeparator("\t");
      reader.SetAddCommandLine(false);
      reader.SetInputPath(m_inpath);
      std::string infile(m_infile);
      if (m_infile.find('|')!=std::string::npos)
        infile=infile.substr(0,infile.find('|'));
      reader.SetInputFile(infile+"|BEGIN_SELANA|END_SELANA");
      reader.SetComment("#");
      m_pt = reader.GetValue<double>("PT_GAMMA", 10);
      m_dr_gamma_lep = reader.GetValue<double>("DELTAR_GAMMA_LEPTON", 0.2);
      m_eps = reader.GetValue<double>("ISO_EPSILON", 0.1);
      m_dr = reader.GetValue<double>("ISO_DELTAR", 0.4);
      m_n = reader.GetValue<double>("ISO_EXPONENT", 1.0);
      m_check_cuts = reader.GetValue<int>("CHECKCUTS", 0);

      msg_Debugging()<<METHOD<<"(): { setting cut variables \n" <<
                       "dr(gamma, lepton)= " << m_dr_gamma_lep <<
                       " pT(gamma) = " << m_pt <<
                       " (dr, eps, n) = (" << m_dr << ", " << m_eps
                    << ", " << m_n <<
                       " ) \n }" << std::endl;

      return true;
    }

    bool Run(ATOOLS::Blob_List *const bl){
      m_direct=false;
      //check, if this event has direct photons TODO: do only once?
      ATOOLS::Blob *sp(bl->FindFirst(ATOOLS::btp::Signal_Process));
      ATOOLS::Particle_Vector v_particles = sp->GetOutParticles();
      for (size_t i(0); i < v_particles.size(); i++){
          ATOOLS::Particle *particle(v_particles.at(i));
          if (particle->Flav().IsPhoton())  m_direct =true;
        }
      v_particles.clear();
      //identify all relevant particles which are active
      ATOOLS::Particle_List particles(bl->ExtractParticles(1));
      ATOOLS::Particle_Vector leptonen, partonen;
      ATOOLS::Particle *leadingphoton;
      double gpt=0;
      for (ATOOLS::Particle_List::const_iterator it=particles.begin(); it!=particles.end(); it++){
          //get photon with highest pT, which comes either from ME or from YFS, but not from hadrons
          ATOOLS::Particle *particle=*it;
          if (particle->Flav().IsPhoton() && particle->ProductionBlob()->Type()==ATOOLS::btp::QED_Radiation){
              if (particle->Momentum().PPerp()>gpt){
                  gpt = particle->Momentum().PPerp();
                  leadingphoton = particle;
                }
            }
          //get vector of all leptons, which come not from hadron decays
          if(particle->Flav().IsLepton() && particle->Flav().Charge()!=0 &&
                                            particle->ProductionBlob()->Type()==ATOOLS::btp::QED_Radiation){
              leptonen.push_back(particle);
            }

          //get vector of partonen or hadrons
          if(particle->Flav().IsQuark() || particle->Flav().IsGluon() || particle->Flav().IsHadron()){
              partonen.push_back(particle);
            }
        }
       msg_Debugging() <<  METHOD << "()" <<  "   NEW EVENT   \n" << 
                   "leading photon:  "  <<  *leadingphoton  << "\n" <<
                   "lepton1: "   << *leptonen[0] << "\n" <<
                   "lepton2: "   << *leptonen[1] << "\n" <<  std::endl;
        
       
      if (leadingphoton){
          bool lep_cut = LeptonCut(*leadingphoton, leptonen);
          bool gamma_cut = PhotonPtCut(*leadingphoton);
          bool iso_cut = IsoCut(*leadingphoton, partonen);
          //direct case: all cuts have to be fullfilled
          //indirect case: at least one cut is not fullfilled
          // indirect, check cuts: no cut is allowed to be fullfilled, use only for debugging
          if (m_direct && lep_cut && gamma_cut && iso_cut) return true;
          if (!m_check_cuts &&!m_direct &&(!lep_cut || !gamma_cut || !iso_cut)) return true;
          if (m_check_cuts &&!m_direct &&(!lep_cut && !gamma_cut && !iso_cut)) return true;

        }
      else return true;

      return false;
    }

    bool Finish(){}

  };// end of class SELANA



} // end of namespace SELAN



using namespace SHERPA;
using namespace SELAN;

DECLARE_GETTER(SELANA,"SELANA",
               Analysis_Interface,Analysis_Arguments);

Analysis_Interface *ATOOLS::Getter
<Analysis_Interface,Analysis_Arguments,SELANA>::
operator()(const Analysis_Arguments &args) const
{
  msg_Info()<<METHOD<<"(): {}"<<std::endl;
  return new SELANA(args.m_inpath,args.m_infile,args.m_outpath);
}

void ATOOLS::Getter<Analysis_Interface,Analysis_Arguments,
SELANA>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Selector Analysis";
}

