#include "MonopoleOutput.h"
#include "crpropa/Units.h"
#include "crpropa/Version.h"
#include "crpropa/Random.h"
#include "crpropa/base64.h"

#include "kiss/string.h"

#include "kiss/string.h"

#include <cstdio>
#include <stdexcept>
#include <iostream>

#ifdef CRPROPA_HAVE_ZLIB
#include <izstream.hpp>
#include <ozstream.hpp>
#endif

namespace crpropa {

// MonopoleOutput --------------------------------------------------------
MonopoleOutput::MonopoleOutput(const std::string &filename) : outfile(filename.c_str(),
				std::ios::binary), out(&outfile), filename(
				filename), count(0), idx(0), lengthScale(kpc), energyScale(EeV) {
	if (!outfile.is_open())
		throw std::runtime_error(std::string("Cannot create file: ") + filename);
}


void MonopoleOutput::printHeader() const {
	*out << "#";
	//Idx Column
	//*out << "\tI";
	//Serial Number Column
	*out << "\tSN";
	//ID Column
	*out << "\tID";
	//Mass column
	*out << "\tM";
	//Trajectory length column
	*out << "\tD";
	//energy
	*out << "\tE";
	//Position
	*out << "\tX\tY\tZ";
	//Direction
	*out << "\tPx\tPy\tPz";
	//Source energy
	*out << "\tE0";
	//Source position
	*out << "\tX0\tY0\tZ0";
	//Source direction
	*out << "\tPx0\tPy0\tPz0";
	*out << "\n";
}

void MonopoleOutput::setEnergyScale(double scale) {
	energyScale = scale;
}

void MonopoleOutput::setLengthScale(double scale) {
	lengthScale = scale;
}

void MonopoleOutput::process(Candidate *c) const {
	char buffer[1024];
	size_t p = 0;
	
	std::locale old_locale = std::locale::global(std::locale::classic());

	/*//idx
	p += std::sprintf(buffer + p, "%i\t",
			idx);*/
			
	//serial number
	p += std::sprintf(buffer + p, "%10lu\t",
			c->getSerialNumber());
			
	//Current ID
	p += std::sprintf(buffer + p, "%10i\t", c->current.getId());		
	
	//Mass
	MCandidate *Mc = MCandidate::convertToMCandidate(c);
	p += std::sprintf(buffer + p, "%8.5E\t", Mc->Mcurrent.getMass()/(gigaelectronvolt/c_squared));
	/*
	if (isDyon(c->current.getId())) {
		MCandidate *candidate = MCandidate::convertToMCandidate(c);
		p += std::sprintf(buffer + p, "%8.5E\t", MCandidate->current.pmass/(gigaelectronvolt/c_squared));
	}
	else {
		p += std::sprintf(buffer + p, "%10i\t", c->current.getMass()/(gigaelectronvolt/c_squared));
	}
	*/
	
	//trajectory length
	p += std::sprintf(buffer + p, "%8.5E\t",
			c->getTrajectoryLength() / lengthScale);
	
	//Energy	
	p += std::sprintf(buffer + p, "%8.5E\t",
		c->current.getEnergy() / energyScale);
	
	//Position
	const Vector3d pos = c->current.getPosition() / lengthScale;
	p += std::sprintf(buffer + p, "%8.5E\t%8.5E\t%8.5E\t", pos.x, pos.y,
			pos.z);
	
	//Direction				
	const Vector3d dir = c->current.getDirection();
	p += std::sprintf(buffer + p, "%8.5E\t%8.5E\t%8.5E\t", dir.x, dir.y,
			dir.z);
			
	//Source energy
	p += std::sprintf(buffer + p, "%8.5E\t",
			c->source.getEnergy() / energyScale);
	
	//Source position
	const Vector3d pos0 = c->source.getPosition() / lengthScale;
	p += std::sprintf(buffer + p, "%8.5E\t%8.5E\t%8.5E\t", pos0.x, pos0.y,
			pos0.z);

	//Source direction
	const Vector3d dir0 = c->source.getDirection();
	p += std::sprintf(buffer + p, "%8.5E\t%8.5E\t%8.5E\t", dir0.x, dir0.y,
			dir0.z);
			
	buffer[p - 1] = '\n';
	
	std::locale::global(old_locale);
	
#pragma omp critical
	{
		if (count == 0)
			printHeader();
		count++;
		out->write(buffer, p);
	}
	
	/*if (!c->isActive()) {
		idx++;
	}*/
}

void MonopoleOutput::close() {
#ifdef CRPROPA_HAVE_ZLIB
	zstream::ogzstream *zs = dynamic_cast<zstream::ogzstream *>(out);
	if (zs) {
		zs->close();
		delete out;
		out = 0;
	}
#endif
	outfile.flush();
}

MonopoleOutput::~MonopoleOutput() {
	close();
}

// ObserverEnergy --------------------------------------------------------
ObserverEnergy::ObserverEnergy(double energyThreshold, double distanceThreshold = 1*kpc) : energyThreshold(energyThreshold), distanceThreshold(distanceThreshold) { }

DetectionState ObserverEnergy::checkDetection(Candidate *candidate) const
{
		double currentDistance = candidate -> getTrajectoryLength();
		double previousEnergy = candidate -> previous.getEnergy();
		double currentEnergy =  candidate -> current.getEnergy();

		if ((currentDistance > distanceThreshold) & (currentEnergy < energyThreshold) & (currentEnergy < previousEnergy))
			return DETECTED;
		else
			return NOTHING;
}

void ObserverEnergy::setEnergyThreshold(double newEnergyThreshold) {
	energyThreshold = newEnergyThreshold;
}

void ObserverEnergy::setDistanceThreshold(double newDistanceThreshold) {
	distanceThreshold = newDistanceThreshold;
}

std::string ObserverEnergy::getDescription() const {
	std::stringstream ss;
	ss << "ObserverEnergy: " << energyThreshold / EeV << " EeV, " << distanceThreshold / kpc << " kpc";
	return ss.str();
}

} // namespace crpropa
