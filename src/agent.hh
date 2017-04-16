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

class Agent {
public:
  uint32_t id = 0;
  Agent* partner;
  Agent* infector = NULL;

  double singlePeriodDeviation;
  double relationshipPeriodDeviation;
  double relationshipChangeDate;
  double age;
  double desired_age;
  double weight; // Used by some pair-matching algorithms

  unsigned short sex;
  unsigned short sexual_orientation;
  bool infected;
  bool initial_relationship;

  unsigned num_partners = 0;
  unsigned num_infected = 0;


  void setRelationshipPeriod(const double currentDate,
                             const DblMatrix& shapeRelationshipPeriod,
                             const DblMatrix& scaleRelationshipPeriod,
                             const double scaleModifierRelationshipPeriod)
  {
    double shape, scale;
    size_t index1, index2;

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
    scale += relationshipPeriodDeviation + partner->relationshipPeriodDeviation;
    scale *= scaleModifierRelationshipPeriod;
    std::weibull_distribution<double> dist(shape, std::max(0.0001, scale) );
    double relationshipLength = dist(rng);
    relationshipChangeDate = std::max(currentDate,
                                      currentDate + relationshipLength);
    partner->relationshipChangeDate = relationshipChangeDate;
  }

  void setSinglePeriod(const double currentDate,
                       const double shape,
                       const double scale,
                       const DblMatrix& probZeroSinglePeriod,
                       const double probZeroSinglePeriodScale)
  {
    std::uniform_real_distribution<double> uni;
    unsigned index = std::min( (unsigned) age - MIN_AGE,
                               (unsigned) MAX_AGE_CSV);

    double prob = probZeroSinglePeriod[index][sex] * probZeroSinglePeriodScale;
    if (uni(rng) < prob) {
      relationshipChangeDate = currentDate;
    } else {
      std::weibull_distribution<double> dist(shape,
                                             std::max(scale +
                                                      singlePeriodDeviation,
                                                      0.0001));
      relationshipChangeDate = std::max(currentDate,
                                        currentDate + (double) dist(rng) * DAY);
    }
  }

  bool isMatchable(const double currentDate) const
  {
    if (partner == NULL && (currentDate + DAY / 2.0) >= relationshipChangeDate)
      return true;
    else
      return false;
  }

  void print(FILE *f = stdout) const
  {
    fprintf(f, "ID,%7u,Age,%3.2f,Sex,%c,Orientation,%c,Desired,%.0f",
            id, age, (sex == MALE ? 'M' : 'F'),
            (sexual_orientation == HETEROSEXUAL ? 'S' : 'G'), desired_age);
    //fprintf(f, ",Risk,%.3f,%.3f,Partner",
    //        relationship_length_factor, binomial_p_relationship_wait);
    if (partner)
      fprintf(f,",%7u", partner->id);
    else
      fprintf(f, ",      0");
    fprintf(f, ",Date,%.3f,Infection,%d\n", relationshipChangeDate, infected);
  }
};

typedef std::vector<Agent *> AgentVector;

void printAgents(const AgentVector&, unsigned, double, FILE *);
std::ostream& operator<<(std::ostream& os, const Agent& agent);

class Partnerships {
public:
  void insert(const uint32_t id1, const uint32_t id2)
  {
    partnerships.insert(combine(id1, id2));
  }

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

  //private:
  uint64_t combine(uint32_t id1, uint32_t id2) const
  {
    uint64_t A = id1 < id2 ? id1 : id2;
    uint64_t B = id1 > id2 ? id1 : id2;
    uint64_t C = A << 32 | B;
    return C;
  }
  std::unordered_set<uint64_t> partnerships;
};


struct PartnershipScore {
  std::vector<Agent *>::iterator partner;
  double score;
};

#endif
