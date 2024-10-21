#ifndef MONOPOLE_OUTPUT_H
#define MONOPOLE_OUTPUT_H

#include "crpropa/module/ParticleCollector.h"
#include "Monopole.h"
#include "crpropa/module/Observer.h"

#include <fstream>

namespace crpropa {
/**
 * \addtogroup Output
 * @{
 */

/**
 @class MonopoleOutput
 @brief Plain text output for monopole information intended for galactic trajectories.
 */
class MonopoleOutput: public Module {
private:
	void printHeader() const;
	std::ostream *out;
	std::ofstream outfile;
	std::string filename;
	mutable size_t count;
	mutable size_t idx;
	double lengthScale;
	double energyScale;

public:
	/** Constructor
	*/
	MonopoleOutput(const std::string &filename);
	~MonopoleOutput();
	/** Set energy scale.
	 @param scale	energy scale (scale = 1 corresponds to 1 Joule)
	 */
	void setEnergyScale(double scale);
	/** Set length scale.
	 @param scale	length scale (scale = 1 corresponds to 1 meter)
	 */
	void setLengthScale(double scale);
	
	void process(Candidate *candidate) const;
	
	void close();
};
/** @}*/


/**
 @class ObserverEnergy
 @brief Detects particles below a certain energy
 */
class ObserverEnergy: public ObserverFeature {
private:
	double energyThreshold;
	double distanceThreshold;
public:
	/** Constructor
	 @param energyTheshold
	*/
	ObserverEnergy(double energyThreshold, double distanceThreshold);
	DetectionState checkDetection(Candidate *candidate) const;
	void setEnergyThreshold(double newEnergyThreshold);
	void setDistanceThreshold(double newDistanceThreshold);
	std::string getDescription() const;
};

} // namespace crpropa

#endif // MONOPOLE_OUTPUT_H
