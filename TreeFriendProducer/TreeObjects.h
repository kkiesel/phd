#include <TLorentzVector.h>
#include <TVector3.h>

#ifndef TreeObjects_h
#define TreeObjects_h

namespace tree {
// In this namespace classes for the trees are defined.

enum genParticles {
	kGenPhoton,
	kGenElectron,
	kGenJet,
	kNearLepton
};

enum jetMatches {
	kJetId,
	kJetPhoton,
	kJetAllPhoton,
	kJetCount
};

enum electronWorkingPoints {
	kNoElectron,
	kVetoElectron,
	kLooseElectron,
	kMediumElectron,
	kTightElectron
};

class Particle {
	public:
		// functions
		float DeltaR( const Particle &p2 ) const;
		float DeltaR( const TLorentzVector &vec2 ) const;
		float DeltaPhi( float phi2 ) const;
		void setStatus( int status );
		bool isStatus( int status ) const;

		// variables
		float pt, eta, phi;
		short bitFlag;
};

class Photon : public Particle {
	public:
		float ptJet() const;
		float _ptJet;
		float _etaJet;
		float _phiJet;
		float sigmaIphiIphi;
		float r9, sigmaIetaIeta, hadTowOverEm;
		float chargedIso, neutralIso, photonIso;
		bool conversionSafeVeto;
		int pixelseed;
		short matchedJetIndex;
};

class Jet : public Particle{
	public:
		//float DeltaR( const Photon &p2 ) const;
		//float DeltaR( const Jet &p2 ) const;
		float bCSV;
		float chargedHadronEnergy,
			neutralHadronEnergy,
			photonEnergy,
			electronEnergy,
			muonEnergy,
			HFHadronEnergy,
			HFEMEnergy,
			chargedEmEnergy,
			chargedMuEnergy,
			neutralEmEnergy;
};

bool EtGreater(const tree::Particle, const tree::Particle);

} // end namespace definition

#endif



