#include "agent.hh"

void printAgents(const AgentVector& agents,
                 unsigned simulation_num,
                 double currentDate,
                 FILE *out)
{
  for (auto& agent: agents) {
    fprintf(out, "AGENT,%u,%f,", simulation_num, currentDate);
    agent->print(out);
  }
}

// We don't use the ostream operator here, preferring the C stdio library
// because it is implemented much faster by the GNU compiler.

std::ostream& operator<<(std::ostream& os, const Agent& agent)
{
  os << "ID," << agent.id << ",Age," << agent.age
     << ",Sex," << (agent.sex == MALE ? "M" : "F")
     << ",Orientation," << (agent.sexual_orientation == HETEROSEXUAL
                            ? "S" : "G")
     << ",Desired," << agent.desired_age;
  // std::ios_base::fmtflags oldflags = os.flags();
  // std::streamsize oldprecision = os.precision();
  // os << std::fixed << std::setprecision(3)
  //    << ",Risk," << agent.relationship_length_factor
  //    << "," << agent.binomial_p_relationship_wait;
  // os.flags (oldflags);
  // os.precision (oldprecision);
  if (agent.partner) {
    auto &p = agent.partner;
    os << ",ID," << p->id << ",Age," << p->age
       << ",Sex," << (p->sex == MALE ? "M" : "F")
       << ",Orientation," << (p->sexual_orientation == HETEROSEXUAL
                              ? "S" : "G")
       << ",Desired," << p->desired_age;
  }
  os << ",Date," << agent.relationshipChangeDate;
  return os;
}
