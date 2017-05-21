#ifndef AGENT_HH
#define AGENT_HH

#include <cstdio>
#include <cstdint>
#include <cstdlib>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <vector>

#include "common.hh"
#include "linear.hh"

/**
   Manages the behaviour of individual agents.
 */

class Agent {
public:
  uint32_t id = 0;
  Agent* partner;
  Agent* infector = NULL;

  double singlePeriodDeviation;
  double relationshipPeriodDeviation;
  double casualSexDeviation;
  double relationshipChangeDate = 0.0;
  double age;
  double desiredAge;
  double weight; // Used by some pair-matching algorithms

  unsigned short sex;
  unsigned short sexualOrientation;
  bool infected;
  bool virgin = false;
  bool casualSex = false;

  unsigned numPartners = 0;
  unsigned numInfected = 0;


  /**
     Sets the length of the period for which two agents maintain a
     relationship before breaking up.

     @param currentDate[in] Current date in the simulation. This is the earliest
     date that the relation can end.
     @param shapeRelationshipPeriod[in] Weibull shape parameter matrix
     @param scaleRelationshipPeriod[in] Weibull scale parameter matrix
     @param scaleModifierRelationshipPeriod[in] Fitting parameter with which the
     Weibull scale parameter is multiplied.
  */
  void setRelationshipPeriod(const double currentDate,
                             const DblMatrix& shapeRelationshipPeriod,
                             const DblMatrix& scaleRelationshipPeriod,
                             const double scaleModifierRelationshipPeriod)
  {
#ifdef LOGGING
    double shape = 0.0, scale = 0.0;
    size_t index1 = 0, index2 = 0;
#else
    double shape, scale;
    size_t index1, index2;
#endif

    if (casualSex == true || partner->casualSex == true) {
      relationshipChangeDate = currentDate;
      partner->relationshipChangeDate = currentDate;
    } else {
      if (sex == MALE) {
        index1 = std::min( (unsigned) age - MIN_AGE, (unsigned) MAX_AGE_CSV);
        index2 = std::min( (unsigned) partner->age - MIN_AGE,
                           (unsigned) MAX_AGE_CSV);
      } else {
        index2 = std::min( (unsigned) age - MIN_AGE, (unsigned) MAX_AGE_CSV);
        index1 = std::min( (unsigned) partner->age - MIN_AGE,
                           (unsigned) MAX_AGE_CSV);
      }

      shape = shapeRelationshipPeriod[ index1 ] [ index2 ];
      scale = scaleRelationshipPeriod[ index1 ] [ index2 ];
      std::weibull_distribution<double> dist(shape, scale);
      double deviation = (relationshipPeriodDeviation +
                          partner->relationshipPeriodDeviation) / 2.0;
      double relationshipLength = dist(rng) * scaleModifierRelationshipPeriod * deviation;
      relationshipChangeDate = std::max(currentDate, currentDate + relationshipLength);
      partner->relationshipChangeDate = relationshipChangeDate;
      virgin = false;
      partner->virgin = false;
    }

#ifdef LOGGING
    std::ostringstream stream;
    stream << "LOGAGENTREL," << id << "," << age << ","
           << sex << "," << sexualOrientation
           << "," << partner->id << "," << index1 << "," << index2 << ","
           << shape << "," << scale  << ","
           << currentDate << "," << relationshipChangeDate << std::endl;
    stream << "LOGAGENTREL," << partner->id << "," << partner->age << ","
           << partner->sex << "," << partner->sexualOrientation
           << "," << id << "," << index1 << "," << index2 << ","
           << shape << "," << scale  << ","
           << currentDate << "," << relationshipChangeDate << std::endl;
    std::cout << stream.str();
#endif

  }

  /**
     Sets the period that an agent stays single before forming a relationship
     again.

     @param currentDate[in] Current date in the simulation. This is the earliest
     date that the relation can end.
     @param weibullSinglePeriod[in] Weibull shape and scale parameters matrix
     @param scaleModifierSinglePeriod[in] Fitting parameter with which the
     Weibull scale parameter is multiplied.
     @param probZeroSinglePeriod[in] Probability agent remains single for zero
     days.
     @param probZeroSinglePeriodScale[in] Fitting parameter with which the
     probability of zero period parameter is multiplied.
   */

  void setSinglePeriod(const double currentDate,
                       const DblMatrix& weibullSinglePeriodFirstTime,
                       const double& scaleVirginPeriod,
                       const DblMatrix& weibullSinglePeriodSubsequentTimes,
                       const double scaleSinglePeriod,
                       const DblMatrix& probZeroSinglePeriod,
                       const double probZeroSinglePeriodScale)
  {
    unsigned index1;
    double shape = 0.0, scale = 0.0;
    std::uniform_real_distribution<double> uni;

    // Process virgins first
    if (relationshipChangeDate == 0.0) {
      index1 = std::min( (size_t) age - MIN_AGE,
                         (size_t) weibullSinglePeriodFirstTime.size() - 1);
      shape = weibullSinglePeriodFirstTime[index1][sex * 2];
      scale = weibullSinglePeriodFirstTime[index1][sex * 2 + 1];
      std::weibull_distribution<double> dist(shape, scale);
      double ageFirstSex  = dist(rng) * scaleVirginPeriod;
      if (ageFirstSex > age) {
        virgin = true;
        double timeFirstSex = ageFirstSex - age;
        relationshipChangeDate = currentDate + timeFirstSex;
      }
    }

    // Process agents who have had sex
    if (virgin == false) {
      index1 = std::min( (unsigned) age - MIN_AGE,
                         (unsigned) MAX_AGE_CSV);
      // Go immediately into a new relationship
      double prob = probZeroSinglePeriod[index1][sex] * probZeroSinglePeriodScale;
      if (uni(rng) < prob) {
        relationshipChangeDate = currentDate;
      } else {
        shape = weibullSinglePeriodSubsequentTimes[index1][sex * 2];
        scale = weibullSinglePeriodSubsequentTimes[index1][sex * 2 + 1];
        std::weibull_distribution<double> dist(shape, scale);
        double timeSingle  = dist(rng) * scaleSinglePeriod * singlePeriodDeviation;
        relationshipChangeDate = currentDate + timeSingle;
      }
    }
#ifdef LOGGING
    std::ostringstream stream;
    stream << "LOGAGENTSIN," << id << "," << age
           << "," << sex << "," << sexualOrientation
           << "," << 0 << "," << index1 << ","
           << shape << "," << scale  << ","
           << currentDate << "," << relationshipChangeDate << std::endl;
    std::cout << stream.str();
#endif
  }

  /**
     Checks if the agent is to be placed in the mating pool.

     @param currentDate[in] Current date of the simulation.
   */

  bool isMatchable(const double currentDate, const DblMatrix& probCasualSex)
  {
    if (partner == NULL) {
      if ( (currentDate + DAY / 2.0) >= relationshipChangeDate) {
        casualSex = false;
        return true;
      } else {
        if (virgin == false) {
          std::uniform_real_distribution<double> uni(0.0, 1.0);
          size_t index = std::min( (size_t) age - MIN_AGE,
                                     (size_t) probCasualSex.size() - 1);
          if (uni(rng) < (probCasualSex[index][1] * casualSexDeviation)) {
            casualSex = true;
            return true;
          }
        }
      }
    }
    return false;
  }
};

typedef std::vector<Agent *> AgentVector;
void printAgents(const AgentVector&, unsigned, double,
                 std::ostream& = std::cout);
std::ostream& operator<<(std::ostream& os, const Agent& agent);

/**
   Tracks partnerships across the simulation. This class is optimised to use
   as little space as possiblee. It simply keeps the IDs of the two partner
   agents in an unordered set (i.e. a hash table) where each element is
   a 64 bit word.
 */

class Partnerships {
public:

  /**
     Inserts a new partnership into the set.

     @param id1[in] ID of one of the agents in the partnership
     @param id2[in] ID of the other agent in the partnership
   */
  void insert(const uint32_t id1, const uint32_t id2)
  {
    partnerships.insert(combine(id1, id2));
  }


  /**
     Checks if a partnership exists.

     @param id1[in] ID of one of the agents in the partnership
     @param id2[in] ID of the other agent in the partnership
  */
  bool exists(const uint32_t id1, const uint32_t id2) const
  {
    if (partnerships.find(combine(id1, id2)) == partnerships.end())
      return false;
    else
      return true;
  }
  size_t size() const
  {
    return partnerships.size();
  }

  /**
     Combine two 32-bit IDs into one 64-bit word for insertion into the hash
     table.

     @param id1[in] ID of one of the agents in the partnership
     @param id2[in] ID of the other agent in the partnership
  */
  uint64_t combine(uint32_t id1, uint32_t id2) const
  {
    uint64_t A = id1 < id2 ? id1 : id2;
    uint64_t B = id1 > id2 ? id1 : id2;
    uint64_t C = A << 32 | B;
    return C;
  }

  /**
     Hash table that stores the partnerships
  */
  std::unordered_set<uint64_t> partnerships;
};


struct PartnershipScore {
  std::vector<Agent *>::iterator partner;
  double score;
};

#endif
