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
  unsigned short sex;
  unsigned short sexualOrientation;
  Agent* partner;
  double age;
  double desiredAge;
  double weight; // Used by some pair-matching algorithms

  Agent* infector = NULL;
  double singlePeriodFactor;
  double relationshipPeriodFactor;
  double casualSexFactor;
  double relationshipChangeDate = 0.0;


  bool infected;
  bool virgin = false;
  bool casualSex = false;

  unsigned numPartners = 0;
  unsigned numInfected = 0;

};

typedef std::vector<Agent *> AgentVector;
void printAgents(const AgentVector&, std::string, unsigned, double,
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
