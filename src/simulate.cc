/**
   stisimulator: simulate.cc
   Purpose: Microsimulation of sexually transmitted infection epidemics

   @author Nathan Geffen
   @version 0.1 23/2/2017
   @license GPL v3
*/

#include <unistd.h>
#include "simulate.hh"

thread_local RNG_TYPE rng;

/**
   Writes analytical information in comma separated format to stdout. This is
   the final function called by other csvout functions.

   @param desc1[in] Description 1
   @param desc2[in] Description 2
   @param value[in] Parameter to print
*/

inline void writeCsvLine(const std::string& simulationName,
                         const unsigned simulationNum,
                         const double date,
                         const std::string& desc1,
                         const std::string& desc2,
                         const std::ostringstream& extra)
{
  // Must do it this way to ensure multithreading doesn't
  // result in one line's output spread over multiple lines
  std::ostringstream stream;
  stream << simulationName << "," << simulationNum << ","
         << std::fixed << std::setprecision(3) << date << ","
         << desc1 << "," << desc2 << ","
         << extra.str() << std::endl;
  std::cout << stream.str();
}


/**
   Calculates the number of people who are single in the initial population.

   @param data[in] matrix of probabilities
   @return number of people who are single
*/

double calcNumberSingles(const DblMatrix& data, const unsigned X)
{
  auto ageshare = getCol(data, 1);
  auto femratio = getCol(data, 2);
  auto relm_share = getCol(data, 3);
  auto relw_share = getCol(data, 4);

  auto t1 = subVector(1, relw_share);
  auto t2 = multVectors(femratio, t1);
  auto t3 = multVectors(t2, ageshare);

  auto t4 = subVector(1, femratio);;
  auto t5 = subVector(1, relm_share);
  auto t6 = multVectors(t4, t5);
  auto t7 = multVectors(t6, ageshare);

  auto t8 = addVector(t3, t7);
  auto t9 = sumVector(t8);

  double num_agents = X * t9;
  return round(num_agents);
}

// Initialization routines

void setInitialInfection(Agent &agent,
                         const std::vector<double>& initialInfectionRatesMSW,
                         const std::vector<double>& initialInfectionRatesMSM,
                         const std::vector<double>& initialInfectionRatesWSM,
                         const std::vector<double>& initialInfectionRatesWSW)
{
  std::uniform_real_distribution<double> uni;
  if (agent.sex == MALE && agent.sexualOrientation == HETEROSEXUAL) {
    agent.infected = (uni(rng) < initialInfectionRatesMSW[(int) agent.age])
      ? true : false;
  } else if (agent.sex == MALE && agent.sexualOrientation == HOMOSEXUAL) {
    agent.infected = (uni(rng) < initialInfectionRatesMSM[(int) agent.age])
      ? true : false;
  } else if (agent.sex == FEMALE && agent.sexualOrientation == HETEROSEXUAL) {
    agent.infected = (uni(rng) < initialInfectionRatesWSM[(int) agent.age])
      ? true : false;
  } else {
    agent.infected = (uni(rng) < initialInfectionRatesWSW[(int) agent.age])
      ? true : false;
  }
}

/****************************/

void ageEvent(Simulation* simulation)
{
  if ( (simulation->inStabilizationPeriod == false) ) {
    for (auto& agent : simulation->agents) {
      agent->age += simulation->timeStep;
    }
  }
}

void infectEvent(Simulation* simulation)
{
  if ( simulation->inStabilizationPeriod == false) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (auto& agent: simulation->agents) {
      if (agent->partner && agent->infected == false &&
          agent->partner->infected == true) {
        double risk_infection;

        if (agent->partner->sex == agent->sex) {
          if (agent->sex == MALE) {
            risk_infection = simulation->homMaleInfectiousness;
          } else {
            risk_infection = simulation->homFemaleInfectiousness;
          }
        } else {
          if (agent->sex == MALE) {
            risk_infection = simulation->hetFemaleInfectiousness;
          } else {
            risk_infection = simulation->hetMaleInfectiousness;
          }
        }

        if (dist(rng) < risk_infection) {
          agent->infected = true;
          agent->infector = agent->partner;
          ++agent->partner->numInfected;
          simulation->trackRiskFactors(agent);
        }
      }
    }
  }
}

void breakupCommon(Simulation *simulation, unsigned breakups)
{
  if (simulation->printNumBreakups) simulation->csvout("BREAKUPS", "", breakups);
  simulation->totalBreakups += breakups;
}

void frequencyBreakupEvent(Simulation* simulation)
{
  unsigned breakups = 0;
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  for (auto& agent: simulation->agents) {
    if (agent->partner && agent->id < agent->partner->id) {
      bool executeBreakup = false;
      if (agent->casualSex == true || agent->partner->casualSex == true) {
        agent->casualSex = agent->partner->casualSex = false;
        ++simulation->casualSexEncounters;
        executeBreakup = true;
      } else {
        size_t index = std::min( (size_t) agent->age - MIN_AGE,
                                 (size_t) simulation->breakupProb.size() - 1);
        double prob1 = simulation->breakupProb[index][1 + agent->sex] *
          agent->relationshipPeriodFactor;
        index = std::min( (size_t) agent->partner->age - MIN_AGE,
                          (size_t) simulation->breakupProb.size() - 1);
        double prob2 = simulation->breakupProb[index][1 + agent->partner->sex] *
          agent->partner->relationshipPeriodFactor;
        double prob = (prob1 + prob2) / 2.0;
        if (uni(rng) < prob) executeBreakup = true;
      }
      if (executeBreakup == true) {
        agent->partner->partner = NULL;
        agent->partner = NULL;
        ++breakups;
      }
    }
  }
  breakupCommon(simulation, breakups);
}

void randomBreakupEvent(Simulation *simulation)
{
  unsigned breakups = 0;
  double mean = simulation->meanRatePairsTimeStep * simulation->agents.size();
  double sd = simulation->sdRatePairs * mean;
  std::normal_distribution<double> norm(mean, sd);
  unsigned maxBreakups = (unsigned) std::max(0.0, norm(rng));
  std::shuffle(simulation->agents.begin(), simulation->agents.end(), rng);
  std::uniform_real_distribution<double> uni(0.0, 1.0);

  auto& agents = simulation->agents;
  for (auto it = agents.begin();
       it != agents.end() && breakups < maxBreakups; ++it) {
    auto agent = *it;
    if (agent->partner) {
      double relationshipPeriodFactor =
        (agent->relationshipPeriodFactor +
         agent->partner->relationshipPeriodFactor) / 2.0;
      if (uni(rng) < relationshipPeriodFactor) {
        agent->partner->partner = NULL;
        agent->partner = NULL;
        ++breakups;
      }
    }
  }
  breakupCommon(simulation, breakups);
}

void limitFrequencyBreakupEvent(Simulation *simulation)
{
  AgentVector partners;
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  // Find all agents who are in partnerships and, to prevent duplication,
  // have a lower id than their partner
  for (auto& agent: simulation->agents) {
    if (agent->partner && agent->id < agent->partner->id) {
      if (agent->casualSex == true) { // Can assume partner always has same value for this
        agent->weight = -1;
      } else {
        size_t index = std::min( (size_t) agent->age - MIN_AGE,
                                 (size_t) simulation->breakupProb.size() - 1);
        double prob1 = simulation->breakupProb[index][1 + agent->sex];
        index = std::min( (size_t) agent->partner->age - MIN_AGE,
                          (size_t) simulation->breakupProb.size() - 1);
        double prob2 = simulation->breakupProb[index][1 + agent->partner->sex];
        double prob = (prob1 + prob2) / 2.0;
        agent->weight = uni(rng) - prob;
      }
      partners.push_back(agent);
    }
  }
  double mean = simulation->meanRatePairsTimeStep * simulation->agents.size();
  double sd = simulation->sdRatePairs * mean;
  std::normal_distribution<double> norm(mean, sd);
  size_t breakups = std::min( (size_t) std::max(0.0, norm(rng)),
                              partners.size());
  nth_element(partners.begin(), partners.begin() + breakups, partners.end(),
              [](const Agent *a, const Agent *b)
              {
                return a->weight < b->weight;
              });

  for (auto it = partners.begin(); it < partners.begin() + breakups; ++it) {
    auto agent = *it;
    if (agent->casualSex) {
      agent->casualSex = false;
      agent->partner->casualSex = false;
      ++simulation->casualSexEncounters;
    }
    agent->partner->partner = NULL;
    agent->partner = NULL;
  }
  breakupCommon(simulation, breakups);
}


void matingPoolCommon(Simulation* simulation)
{
  if (simulation->matingPool.size() % 2 == 1) simulation->matingPool.pop_back();
  if (simulation->printNumMatings) {
    simulation->csvout("MATINGPOOL", "", simulation->matingPool.size());
  }
}


void frequencyMatingPoolEvent(Simulation* simulation)
{
  simulation->matingPool.clear();
  for (auto& agent: simulation->agents) {
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    if (agent->partner == NULL) {
      size_t index =
        std::min( (size_t) agent->age - MIN_AGE,
                  (size_t) simulation->relationshipProb.size() - 1);
      if (uni(rng) < simulation->relationshipProb[index][1 + agent->sex] *
          agent->singlePeriodFactor) {
          simulation->matingPool.push_back(agent);
      } else {
        size_t index = std::min( (size_t) agent->age - MIN_AGE,
                                 (size_t) simulation->casualSexProb.size() - 1);
        if (uni(rng) < (simulation->casualSexProb[index][1 + agent->sex] *
                        agent->casualSexFactor)) {
          agent->casualSex = true;
          simulation->matingPool.push_back(agent);
        }
      }
      // size_t index = std::min( (size_t) agent->age - MIN_AGE,
      //                          (size_t) simulation->casualSexProb.size() - 1);
      // if (uni(rng) < (simulation->casualSexProb[index][1 + agent->sex] *
      //                 agent->casualSexFactor)) {
      //   agent->casualSex = true;
      //   simulation->matingPool.push_back(agent);
      // } else {
      //   size_t index =
      //     std::min( (size_t) agent->age - MIN_AGE,
      //               (size_t) simulation->relationshipProb.size() - 1);
      //   if (uni(rng) < simulation->relationshipProb[index][1 + agent->sex] *
      //       agent->singlePeriodFactor) {
      //     simulation->matingPool.push_back(agent);
      //   }
      // }

    }
  }
  std::shuffle(simulation->matingPool.begin(), simulation->matingPool.end(),
               rng);
  matingPoolCommon(simulation);
}


void limitFrequencyMatingPoolEvent(Simulation *simulation)
{
  simulation->matingPool.clear();
  std::uniform_real_distribution<double> uni;

  for (auto& agent: simulation->agents) {
    if (agent->partner == NULL) {
      size_t index = std::min( (size_t) agent->age - MIN_AGE,
                               (size_t) simulation->casualSexProb.size() - 1);
      if (uni(rng) < (simulation->casualSexProb[index][1 + agent->sex] *
                      agent->casualSexFactor)) {
        agent->casualSex = true;
        agent->weight = -1.0;
      } else {
        size_t index = std::min( (size_t) agent->age - MIN_AGE,
                                 (size_t) simulation->breakupProb.size() - 1);
        double prob = simulation->relationshipProb[index][1 + agent->sex];
        agent->weight = uni(rng) - prob;
      }
      simulation->matingPool.push_back(agent);
    }
  }

  double mean = simulation->meanRatePairsTimeStep * simulation->agents.size();
  double sd = simulation->sdRatePairs * mean;
  std::normal_distribution<double> norm(mean, sd);
  size_t maxMatingAgents = std::min( (size_t) std::max(0.0, norm(rng)) * 2,
                                     simulation->matingPool.size());

  nth_element(simulation->matingPool.begin(),
              simulation->matingPool.begin() + maxMatingAgents,
              simulation->matingPool.end(),
              [](const Agent *a, const Agent *b)
              {
                return a->weight < b->weight;
              });

  simulation->matingPool.resize(maxMatingAgents);
  std::shuffle(simulation->matingPool.begin(), simulation->matingPool.end(), rng);
  matingPoolCommon(simulation);
}

void randomMatingPoolEvent(Simulation* simulation)
{
  unsigned matingAgents = 0;
  double mean = simulation->meanRatePairsTimeStep * simulation->agents.size();
  double sd = simulation->sdRatePairs * mean;
  std::normal_distribution<double> norm(mean, sd);
  unsigned maxMatingAgents = (unsigned) std::max(0.0, norm(rng)) * 2;
  std::shuffle(simulation->agents.begin(), simulation->agents.end(), rng);

  simulation->matingPool.clear();
  std::uniform_real_distribution<double> uni;
  for (auto it = simulation->agents.begin();
       it != simulation->agents.end() && matingAgents < maxMatingAgents; ++it) {
    auto agent = *it;
    if (agent->partner == NULL && uni(rng) < agent->singlePeriodFactor) {
      simulation->matingPool.push_back(agent);
      ++matingAgents;
    }
  }
  matingPoolCommon(simulation);
}

void randomMatchEvent(Simulation* simulation)
{
  AgentVector& unmatchedAgents = simulation->matingPool;

  if (unmatchedAgents.size() > 0)
    for (size_t i = 0; i < unmatchedAgents.size()  - 1; i += 2)
      simulation->makePartner(unmatchedAgents[i], unmatchedAgents[i + 1],
                              simulation->distance(unmatchedAgents[i],
                                                   unmatchedAgents[i + 1]));
}

void randomKMatchEvent(Simulation *simulation)
{
  AgentVector& unmatchedAgents = simulation->matingPool;

  if(unmatchedAgents.size()) {
    for (auto it = unmatchedAgents.begin(); it < unmatchedAgents.end() - 1;
         ++it) {
      if ( (*it)->partner == NULL) {
        auto last =  (unmatchedAgents.end() - it) < (simulation->neighbors + 1)
          ? unmatchedAgents.end() : it + simulation->neighbors + 1;
        auto partnershipScore = simulation->closestPairMatchN(it, last);
        if (partnershipScore.partner != last)
          simulation->makePartner(*it, *partnershipScore.partner,
                                  partnershipScore.score);
      }
    }
  }
}

void clusterShuffleMatchEvent(Simulation* simulation)
{
  AgentVector& unmatchedAgents = simulation->matingPool;
  uint64_t cluster_size = unmatchedAgents.size() / simulation->clusters;
  for (auto& a : unmatchedAgents) a->weight = simulation->clusterValue(a);
  sort(unmatchedAgents.rbegin(), unmatchedAgents.rend(), [](Agent *a, Agent *b)
       { return a->weight < b->weight; });
  for (uint64_t i = 0; i < simulation->clusters; ++i) {
    auto first = unmatchedAgents.begin() + i * cluster_size;
    auto last = first + cluster_size;
    if (last > unmatchedAgents.end()) last = unmatchedAgents.end();
    std::shuffle(first, last, rng);
  }
  if(unmatchedAgents.size()) {
    for (auto it = unmatchedAgents.begin(); it < unmatchedAgents.end() - 1;
         ++it) {
      if ( (*it)->partner == NULL) {
        auto last = (unmatchedAgents.end() - it) < (simulation->neighbors + 1)
          ? unmatchedAgents.end() : it + simulation->neighbors + 1;
        auto partnershipScore = simulation->closestPairMatchN(it, last);
        if (partnershipScore.partner != last)
          simulation->makePartner(*it, *partnershipScore.partner,
                                  partnershipScore.score);
      }
    }
  }
}


/******************** Call Blossom V reference algorithm ******************/

void graphPairs(const char *graph,
                const Simulation* simulation,
                AgentVector& agents)
{
  FILE *f = fopen(graph, "w");

  uint64_t vertices = (uint64_t) agents.size();
  uint64_t edges = vertices * (vertices - 1) / 2;
  fprintf(f, "%lu %lu\n", vertices, edges);
  for (uint64_t i = 0; i < agents.size(); ++i) {
    for (uint64_t j = i + 1; j < agents.size(); ++j) {
      double d = simulation->distance(agents[i],agents[j]);
      fprintf(f, "%lu %lu %.0f\n", i, j, d * GRAPH_ACCURACY);
    }
  }
  fclose(f);
}

void blossomVMatchEvent(Simulation* simulation)
{
  AgentVector& agents = simulation->matingPool;

  if (agents.size() == 0) return;

  std::ostringstream ss;
  ss << std::this_thread::get_id();
  std::string suffix = ss.str() + std::string(".txt");
  std::string graph_file = std::string("bv_graph_") + suffix;
  std::string blossom_out_file = std::string("bv_out_") + suffix;

  graphPairs(graph_file.c_str(), simulation, agents);
  std::string command("./blossom5 -e ");
  command += graph_file;
  command +=  std::string(" -w ");
  command += blossom_out_file;
  command += std::string(" > /dev/null");
  fflush(stdout);
  if(system(command.c_str()) != 0) {
    throw std::runtime_error("Error executing Blossom V");
  }

  FILE *f = fopen(blossom_out_file.c_str(), "r");
  uint64_t from, to;
  if (fscanf(f, "%lu %lu\n", &from, &to) != 2) {
    fclose(f);
    throw std::runtime_error("Unexpected format in Blossom V graph file");
  }
  for (size_t i = 0; i < agents.size() / 2; ++i) {
    if (fscanf(f, "%lu %lu\n", &from, &to) != 2) {
      throw std::runtime_error("Unexpected format in Blossom V graph file");
    }
    simulation->makePartner(agents[from], agents[to],
                            simulation->distance(agents[from], agents[to]));
  }
  fclose(f);
}

/*********************************************************************/

void Simulation::setEvents()
{
  if (parameterMap.at("AGE_EVENT").isSet()) events.push_back(ageEvent);
  if (parameterMap.at("INFECT_EVENT").isSet()) events.push_back(infectEvent);

  string s = parameterMap.at("BREAKUP_EVENT").str();
  if (s == "FREQUENCY") {
    events.push_back(frequencyBreakupEvent);
  } else if (s == "RANDOM") {
    events.push_back(randomBreakupEvent);
  } else if (s == "LIMIT") {
    events.push_back(limitFrequencyBreakupEvent);
  } else if (s != "NONE") {
    throw std::runtime_error("Unknown breakup event");
  };

  s = parameterMap.at("MATING_POOL_EVENT").str();
  if (s == "FREQUENCY") {
    events.push_back(frequencyMatingPoolEvent);
  } else if (s == "RANDOM") {
    events.push_back(randomMatingPoolEvent);
  } else if (s == "LIMIT") {
    events.push_back(limitFrequencyMatingPoolEvent);
  } else if (s != "NONE") {
    throw std::runtime_error("Unknown breakup event");
  };

  s = parameterMap.at("MATCH_EVENT").str();
  if (s == "RPM") {
    events.push_back(randomMatchEvent);
  } else if (s == "RKPM") {
    events.push_back(randomKMatchEvent);
  } else if (s == "CSPM") {
    events.push_back(clusterShuffleMatchEvent);
  }  else if (s == "BLOSSOMV") {
    events.push_back(blossomVMatchEvent);
  } else {
    throw std::runtime_error("Unknown matching algorithm");
  }
}


void callSimulation(ParameterMap parameterMap,
                    const SymbolTable& symbolTable,
                    const unsigned parmPrint,
                    const unsigned simulationFrom,
                    const unsigned simulationTo)
{

  for (unsigned i = simulationFrom; i < simulationTo; ++i) {
    for (auto& symbol : symbolTable.table[i]) {
      parameterMap[symbol.first].values = {symbol.second};
    }

    std::string prefix;
    std::string suffix = ",";
    prefix.append(parameterMap.at("SIMULATION_NAME").str());
    prefix.append(",");
    // Parameter printing
    if (parmPrint > NO_PARMS) {
      double date = parameterMap.at("START_DATE").dbl();
      std::string name(parameterMap.at("SIMULATION_NAME").str());
      for (auto& p: parameterMap) {
        if (parmPrint == ALL_PARMS || p.second.isVaryingRange == true) {
          std::ostringstream stream;
          if (p.second.isString) {
            stream << p.second.str();
          } else {
            for (size_t i = 0; i < p.second.values.size(); ++i) {
              stream << p.second.values[i];
              if (i + 1 < p.second.values.size()) stream << ",";
            }
          }
          writeCsvLine(name, i, date, "PARAMETER", p.first, stream);
        }
      }
    }
    // Run the simulation
    Simulation(parameterMap, i).simulate();
  }
}


void execSimulationSet(std::vector<ParameterMap> parameterMaps,
                       unsigned seed)
{
  bool csvHeaderWritten = false;

  for (auto& parameterMap: parameterMaps) {
    if (seed == 0) { // Use time
      parameterMap["RANDOM_SEED"].values = { (double) clock() };
    } else if (seed != 1) {
      parameterMap["RANDOM_SEED"].values = { (double) seed};
    }

    // Construct symbol table
    unsigned numSimulations = parameterMap.at("NUM_SIMULATIONS").dbl();
    SymbolTable symbolTable(numSimulations);
    for (auto& p: parameterMap.varyingParameters) {
      auto parameter = parameterMap.at(p);
      symbolTable.addSymbol(p, parameter.parentParameter, parameter.values);
    }

#ifdef DEBUG
    if (parameterMap.at("PRINT_SYMBOL_TABLE_SIZE_AND_EXIT").isSet()) {
      size_t cycleLength = 0;
      for (auto& e: symbolTable.initialValues) {
          if (e.second > cycleLength) cycleLength = e.second;
      }
      std::cerr << "Symbol table cycle length: " << cycleLength << std::endl;
      exit(1);
    }
#endif

    // Write the CSV header
    if (parameterMap.at("CSV_HEADER").isSet() && csvHeaderWritten == false) {
      std::cout << "Name,Num,Date,Desc1,Desc2,Value" << std::endl;
      csvHeaderWritten = true;
    }

    // Needed for parameter printing
    unsigned parmPrint = (unsigned) parameterMap.at("PRINT_PARAMETERS").dbl();

    // Organise threads
    unsigned numThreads = parameterMap.at("NUM_THREADS").dbl();
    if (numThreads == 0) {
      numThreads = sysconf(_SC_NPROCESSORS_ONLN);
      if (numThreads == 0) numThreads = 1;
    }
    if (numSimulations < numThreads) numThreads = numSimulations;
    unsigned simulationsPerThread = numSimulations / numThreads;
    while (simulationsPerThread * numThreads < numSimulations) {
      ++simulationsPerThread;
    }

    std::vector<std::thread> thread(numThreads);

    for (unsigned i = 0; i < numThreads; ++i) {
      unsigned simulationFrom = i * simulationsPerThread;
      unsigned simulationTo = std::min( (i + 1) * simulationsPerThread, numSimulations);
      thread[i] = std::thread(callSimulation, parameterMap,
                              std::ref(symbolTable), parmPrint,
                              simulationFrom, simulationTo);
    }

    for (unsigned i = 0; i < numThreads; ++i)  thread[i].join();
  }
}


/**
   Runs a suite of tests to check that nothing obvious is broken.

   @param parameterMap parameters with which to run the test suite
*/

void runTests(ParameterMap& parameterMap)
{
  unsigned successes = 0, failures = 0;

  // Partnerships
  Partnerships partnerships;
  partnerships.insert(12995, 271);
  partnerships.insert(3, 23994);
  TESTEQ(partnerships.exists(12995, 271), true, successes, failures);
  TESTEQ(partnerships.exists(271, 12995), true, successes, failures);
  TESTEQ(partnerships.exists(3, 23994), true, successes, failures);
  TESTEQ(partnerships.exists(23994, 3), true, successes, failures);
  TESTEQ(partnerships.exists(30, 23994), false, successes, failures);
  TESTEQ(partnerships.exists(23994, 30), false, successes, failures);
  TESTEQ(partnerships.exists(12995, 272), false, successes, failures);
  TESTEQ(partnerships.size() == 2, true, successes, failures);

  Simulation simulation(parameterMap, 0);
  simulation.initializeAgents();

  unsigned expectedAgents = (unsigned) parameterMap.at("NUM_AGENTS").dbl();
  unsigned actualAgents = simulation.agents.size();
  bool correctAgents = (actualAgents == expectedAgents ||
                        actualAgents == expectedAgents - 1) ? true : false;
  TESTEQ(correctAgents, true, successes, failures);

  simulation.agents[0]->age = 20;
  simulation.agents[0]->sex = MALE;
  simulation.agents[0]->desiredAge = 50;
  simulation.agents[0]->sexualOrientation = HETEROSEXUAL;

  simulation.agents[1]->age = 30;
  simulation.agents[1]->sex = MALE;
  simulation.agents[1]->desiredAge = 25;
  simulation.agents[1]->sexualOrientation = HOMOSEXUAL;

  simulation.agents[2]->age = 40;
  simulation.agents[2]->sex = FEMALE;
  simulation.agents[2]->desiredAge = 30;
  simulation.agents[2]->sexualOrientation = HETEROSEXUAL;

  simulation.agents[3]->age = 25;
  simulation.agents[3]->sex = FEMALE;
  simulation.agents[3]->desiredAge = 25;
  // simulation.agents[3]->relationship_length_factor = 0;
  simulation.agents[3]->sexualOrientation = HOMOSEXUAL;

  simulation.agents[4]->age = 30;
  simulation.agents[4]->sex = FEMALE;
  simulation.agents[4]->desiredAge = 20;
  simulation.agents[4]->sexualOrientation = HOMOSEXUAL;

  simulation.agents[5]->age = 30;
  simulation.agents[5]->sex = MALE;
  simulation.agents[5]->desiredAge = 20;
  simulation.agents[5]->sexualOrientation = HOMOSEXUAL;

  simulation.agents[6]->age = 25;
  simulation.agents[6]->sex = MALE;
  simulation.agents[6]->desiredAge = 25;
  // simulation.agents[6]->relationship_length_factor = 0;
  simulation.agents[6]->sexualOrientation = HETEROSEXUAL;


  double d1 = simulation.heuristicDistance(simulation.agents[0],
                                           simulation.agents[2]);
  TESTEQ(d1, 10.0, successes, failures);

  double d2 = simulation.heuristicDistance(simulation.agents[2],
                                           simulation.agents[0]);
  TESTEQ(d1, d2, successes, failures);

  d1 = simulation.heuristicDistance(simulation.agents[0],
                                    simulation.agents[3]);
  TESTEQ(d1, 65.0, successes, failures);
  d2 = simulation.heuristicDistance(simulation.agents[3],
                                    simulation.agents[0]);
  TESTEQ(d1, d2, successes, failures);

  d1 = simulation.heuristicDistance(simulation.agents[0],
                                    simulation.agents[1]);
  TESTEQ(d1, 62.5, successes, failures);
  d2 = simulation.heuristicDistance(simulation.agents[1],
                                    simulation.agents[0]);
  TESTEQ(d1, d2, successes, failures);
  d1 = simulation.heuristicDistance(simulation.agents[1],
                                    simulation.agents[2]);
  TESTEQ(d1, 57.5, successes, failures);
  d2 = simulation.heuristicDistance(simulation.agents[2],
                                    simulation.agents[1]);
  TESTEQ(d1, d2, successes, failures);

  d1 = simulation.heuristicDistance(simulation.agents[1],
                                    simulation.agents[3]);
  TESTEQ(d1, 52.5, successes, failures);
  d1 = simulation.heuristicDistance(simulation.agents[2],
                                    simulation.agents[3]);
  TESTEQ(d1, 60, successes, failures);

  d1 = simulation.heuristicDistance(simulation.agents[3],
                                    simulation.agents[4]);
  TESTEQ(d1, 5, successes, failures);

  d1 = simulation.heuristicDistance(simulation.agents[1],
                                    simulation.agents[5]);
  TESTEQ(d1, 7.5, successes, failures);

  d1 = simulation.tableDistance(simulation.agents[0],
                                simulation.agents[3]);
  TESTEQ(fabs(d1 - 92.85787107631695) < 0.00000000001, true, successes, failures);
  d2 = simulation.tableDistance(simulation.agents[3],
                                simulation.agents[0]);
  TESTEQ(d1, d2, successes, failures);
  d1 = simulation.tableDistance(simulation.agents[3],
                                simulation.agents[4]);
  TESTEQ(fabs(d1 - 21.9126355725) < 0.00000001, true, successes, failures);
  d2 = simulation.tableDistance(simulation.agents[4],
                                simulation.agents[3]);
  TESTEQ(d1, d2, successes, failures);
  d1 = simulation.tableDistance(simulation.agents[1],
                                simulation.agents[5]);
  TESTEQ(fabs(d1 - 0.644530208167975) < 0.00000001, true, successes, failures);
  d2 = simulation.tableDistance(simulation.agents[5],
                                simulation.agents[1]);
  TESTEQ(d1, d2, successes, failures);

  simulation.agents[6]->partner = simulation.agents[3];
  simulation.agents[3]->partner = simulation.agents[6];

  SymbolTable symbolTable(100);
  symbolTable.addSymbol("var1", "", setRange(0, 10, 2));
  symbolTable.addSymbol("var2", "var1", setRange(100.0, 1000.0, 400.0));

  TESTEQ(symbolTable.size(), 100U, successes, failures);
  TESTEQ(symbolTable[0].size(), 2U, successes, failures);
  TESTEQ(symbolTable[99].size(), 2U, successes, failures);
  TESTEQ(symbolTable[0]["var1"], 0.0, successes, failures);
  TESTEQ(symbolTable[0]["var2"], 100.0, successes, failures);
  TESTEQ(symbolTable[6]["var1"], 2.0, successes, failures);
  TESTEQ(symbolTable[6]["var2"], 500.0, successes, failures);
  TESTEQ(symbolTable[99]["var1"], 8.0, successes, failures);
  TESTEQ(symbolTable[99]["var2"], 500.0, successes, failures);

  {
    ParameterMap p;
    Simulation s(p, 1);
    TESTEQ(s.mswAgeDist.size(), ( (size_t) 89), successes, failures);
    TESTEQ( fabs(s.mswAgeDist[78][88] - 0.0039946161 < EPSILON), true, successes, failures);
    truncateMatrix(s.mswAgeDist, 38, 38);
    TESTEQ(s.mswAgeDist.size(), ( (size_t) 38), successes, failures);
    TESTEQ(s.mswAgeDist[0].size(), ( (size_t) 38), successes, failures);
    TESTEQ( fabs(s.mswAgeDist[37][37] - 0.133143 < EPSILON), true, successes, failures);
  }

  {
    parameterMap["START_DATE"].values = {2017};
    parameterMap["END_DATE"].values = {2018};
    parameterMap["TIME_STEP"].values = {0.01};
    parameterMap["HET_MALE_INFECTIOUSNESS"].values = {0.8};
    parameterMap["HET_FEMALE_INFECTIOUSNESS"].values = {0.4};
    parameterMap["PRINT_PARAMETERS"].values = {0};
    parameterMap["ANALYZE_INIT"].values = {0};
    parameterMap["ANALYZE_DURING_SIM"].values = {0};
    parameterMap["ANALYZE_AFTER"].values = {0};
    parameterMap["OUTPUT_NUM_BREAKUPS"].values = {0};
    parameterMap["OUTPUT_NUM_MATINGPOOL"].values = {0};
    parameterMap["OUTPUT_TIMING_AFTER"].values = {0};
    Simulation s(parameterMap, 1);
    s.initializeAgents();
    TESTEQ(s.numInfectedMales > 0, true, successes, failures);
    TESTEQ(s.numInfectedFemales > 0, true, successes, failures);
    unsigned infected = s.numInfectedMales + s.numInfectedFemales;
    TESTEQ(infected < parameterMap.at("NUM_AGENTS").dbl(), true,
           successes, failures);
    unsigned totalinfected = 0;
    for (auto &a : s.agents) if (a->infected) ++totalinfected;
    TESTEQ(infected, totalinfected, successes, failures);
    s.simulate(false);
    totalinfected = 0;
    for (auto &a : s.agents) if (a->infected) ++totalinfected;
    TESTEQ(s.numInfectedMales + s.numInfectedFemales, totalinfected, successes, failures);
    TESTEQ(s.numInfectedMales + s.numInfectedFemales > infected, true,
           successes, failures);
    infected = s.numInfectedMales + s.numInfectedFemales;
    TESTEQ(infected < parameterMap.at("NUM_AGENTS").dbl(), true,
           successes, failures);
  }
  std::cout << "Successes: " << successes
            << " Failures: " << failures << std::endl;
}
