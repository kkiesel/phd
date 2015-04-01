#include "TreeObjects.h"

// Particle

bool tree::EtGreater(const tree::Particle p1, const tree::Particle p2) {
	return p1.pt > p2.pt;
}

float tree::Particle::DeltaR( const Particle &p2 ) const {
	TVector2 tempVect;
	float dphi = tempVect.Phi_mpi_pi( phi - p2.phi );
	float deta = eta - p2.eta;
	return sqrt( dphi*dphi + deta*deta );
}

float tree::Particle::DeltaR( const TLorentzVector &vec2 ) const {
	TVector3 vec1;
	vec1.SetPtEtaPhi( 1, eta, phi );
	return vec1.DeltaR( vec2.Vect() );
}

float tree::Particle::DeltaPhi( float phi2 ) const {
	TVector2 tempVect;
	return tempVect.Phi_mpi_pi( phi - phi2 );
}

bool tree::Particle::isStatus( int status ) const {
	return bitFlag & 1 << status;
}

void tree::Particle::setStatus( int status ) {
	bitFlag |= 1 << status;
}

// Photon

float tree::Photon::ptJet() const {
	return _ptJet ? _ptJet : pt;
}

