#include "forceManager.h"
#include "utilities.h"

using namespace PSim;
using namespace std;

/**
 * @class NonLinear
 * @author Sawyer Hopkins
 * @date 04/26/16
 * @file force.h
 * @brief NonLinear Potential
 */
class NonLinear : public PSim::IForce
{

private:

		//Variables vital to the force.
		double kT;
		double wellDepth;
		int ljNum;
		double charge;
		double salt;
		double cutOff;
		double cutOffSquared;
		double debyeLength; //k
		double mass; // m
		double radius; // r
		bool output;
		long callCount;
		double sphereScale;
		double vP;
		double v0;
		double* volume;
		double* gradient;

	public:

		/**
		 * @brief Creates an new AO Potential.
		 * @param cfg The address of the configuration file reader.
		 */
		NonLinear(config* cfg);
		/**
		 * @brief Releases the force from memory.
		 */
		~NonLinear();

		/**
		 * @brief Get the force from the AO Potential.
		 * @param index The index particle to calculated the force on.
		 * @param nPart The number of particles in the system.
		 * @param boxSize The size of the system.
		 * @param time The current system time.
		 * @param itemCell The cell containing the index particle.
		 * @param items All particles in the system.
		 */
		void getAcceleration(int index, double* sortedParticles, double* particleForce, vector<tuple<int,int>>* particleHashIndex, vector<tuple<int,int>>* cellStartEnd, systemState* state);
		/**
		 * @brief Flag for a force dependent time.
		 * @return True for time dependent. False otherwise.
		 */
		bool isTimeDependent() { return false; }
		/**
		 * @brief Checks for particle interation between the index particle and all particles in the provided cell.
		 * @param boxSize The size of the system.
		 * @param time The current system time.
		 * @param index The particle to find the force on.
		 * @param itemCell The cell to check for interactions in.
		 */
		type3<double> iterCells(int index, int hash, double* sortedParticles, vector<tuple<int,int>>* cellStartEnd, systemState* state);

		void postRoutine(int index, double* sortedParticles, double* particleForce, vector<tuple<int,int>>* particleHashIndex, vector<tuple<int,int>>* cellStartEnd, systemState* state);

		type3<double> postIterCells(int index, int hash, double* sortedParticles, vector<tuple<int,int>>* cellStartEnd, systemState* state);

		void quench(systemState* state);

};

//Class factories.
extern "C" PSim::IForce* getForce(config* cfg)
{
	return new NonLinear(cfg);
}
