#include "agent.hh"

/**
    Prints information about agents to an output file.

    @param agents[in] Agents to print
    @param simulationNum[in] Unique identifier of this simulation
    @param currentDate[in] Current date of the simulation
    @param out[out] File handle for output, defaults to cout
 */

void printAgents(const AgentVector& agents,
                 const unsigned simulationNum,
                 const double currentDate,
                 std::ostream& out)
{
  for (auto& agent: agents) {
    std::ostringstream stream;
    stream << "AGENT," << simulationNum << "," << currentDate << ","
        << *agent << std::endl;
    out << stream.str();
  }
}

/**
   Ostream << for agents.

   @param os[out] Stream to send to
   @param agent[in] Agent to output
*/
std::ostream& operator<<(std::ostream& os, const Agent& agent)
{
  std::ostringstream stream;
  stream << "ID," << agent.id << ",Age," << agent.age
         << ",Sex," << (agent.sex == MALE ? "M" : "F")
         << ",Orientation," << (agent.sexual_orientation == HETEROSEXUAL
                                ? "S" : "G")
         << ",Desired," << agent.desired_age;

  if (agent.partner) {
    auto &p = agent.partner;
    stream << ",ID," << p->id << ",Age," << p->age
           << ",Sex," << (p->sex == MALE ? "M" : "F")
           << ",Orientation," << (p->sexual_orientation == HETEROSEXUAL
                                  ? "S" : "G")
           << ",Desired," << p->desired_age;
  } else {
    stream << ",,,,,,,,,,";
  }
  stream << ",Date," << agent.relationshipChangeDate;
  os << stream.str();
  return os;
}
