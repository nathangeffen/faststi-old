/**
   stisimulator: simulate.cc
   Purpose: Microsimulation of sexually transmitted infection epidemics

   @author Nathan Geffen
   @version 0.1 23/2/2017
   @license GPL v3
*/

#include "simulate.hh"

thread_local std::mt19937 rng;

/**
   Calculates the number of people who are single in the initial population.

   @param data matrix of probabilities
   @return number of people who are single
*/

double calcNumberSingles(DblMatrix data, unsigned X)
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
  if (agent.sex == MALE && agent.sexual_orientation == HETEROSEXUAL) {
    agent.infected = (uni(rng) < initialInfectionRatesMSW[(int) agent.age])
      ? true : false;
  } else if (agent.sex == MALE && agent.sexual_orientation == HOMOSEXUAL) {
    agent.infected = (uni(rng) < initialInfectionRatesMSM[(int) agent.age])
      ? true : false;
  } else if (agent.sex == FEMALE && agent.sexual_orientation == HETEROSEXUAL) {
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
  for (auto& agent : simulation->agents) agent->age += simulation->timeStep;
}

void infectEvent(Simulation* simulation)
{
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for (auto& agent: simulation->agents) {
    if (agent->partner && agent->infected == false &&
        agent->partner->infected == true) {
      double risk_infection;

      if (agent->partner->sex == agent->sex) {
        if (agent->sex == MALE) {
          risk_infection = simulation->hom_male_infectiousness;
        } else {
          risk_infection = simulation->hom_female_infectiousness;
        }
      } else {
        if (agent->sex == MALE) {
          risk_infection = simulation->het_female_infectiousness;
        } else {
          risk_infection = simulation->het_male_infectiousness;
        }
      }

      if (dist(rng) < risk_infection) {
        agent->infected = true;
        agent->infector = agent->partner;
        ++agent->partner->num_infected;
        simulation->trackRiskFactors(agent);
      }
    }
  }
}

void breakupEvent(Simulation* simulation)
{
  unsigned breakups = 0;
  for (auto& agent: simulation->agents) {
    if (agent->partner &&
        ( (simulation->currentDate + DAY / 2.0) >=
          agent->relationshipChangeDate) ) {
      Agent* partner = agent->partner;
      agent->partner = NULL;
      partner->partner = NULL;
      agent->setSinglePeriod(simulation->currentDate,
                             simulation->shapeSinglePeriodDuring,
                             simulation->scaleSinglePeriodDuring);
      partner->setSinglePeriod(simulation->currentDate,
                               simulation->shapeSinglePeriodDuring,
                               simulation->scaleSinglePeriodDuring);
      ++breakups;
    }
  }
  if (simulation->printNumBreakups) {
    printf("%s,BREAKUPS,%u,%.3f,%u\n", simulation->simulationName.c_str(),
           simulation->simulationNum, simulation->currentDate,
           breakups);
  }
  simulation->totalBreakups += breakups;
}

void randomMatchEvent(Simulation* simulation)
{
  // struct timeval timeBegin, timeEnd;
  // double elapsedTime;
  // gettimeofday(&timeBegin, NULL);

  AgentVector unmatchedAgents = simulation->getShuffledUnmatchedAgents();

  if (unmatchedAgents.size() > 0)
    for (size_t i = 0; i < unmatchedAgents.size()  - 1; i += 2)
      simulation->makePartner(unmatchedAgents[i], unmatchedAgents[i + 1],
                              simulation->distance(unmatchedAgents[i],
                                                   unmatchedAgents[i + 1]));
}

void randomKMatchEvent(Simulation *simulation)
{
  std::uniform_real_distribution<double> uni;
  AgentVector unmatchedAgents = simulation->getShuffledUnmatchedAgents();

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
  AgentVector unmatchedAgents = simulation->getShuffledUnmatchedAgents();
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
  AgentVector agents = simulation->getShuffledUnmatchedAgents();

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
  if (parameterMap.at("BREAKUP_EVENT").isSet()) events.push_back(breakupEvent);

  string s = parameterMap.at("MATCH_EVENT").str();
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
                    const unsigned parmPrint,
                    const unsigned simulationNum)
{
  std::string prefix;
  std::string suffix = ",";
  prefix.append(parameterMap["SIMULATION_NAME"].str());
  prefix.append(",");
  prefix.append(std::to_string(simulationNum));
  prefix.append(",");
  // Parameter printing
  if (parmPrint == ALL_PARMS) {
    parameterMap.print(prefix, suffix);
  } else if (parmPrint == VARYING_PARMS) {
    parameterMap.print(prefix, suffix, false, true);
  }

  Simulation(parameterMap, simulationNum).simulate();
}


void execSimulationSet(std::vector<ParameterMap> parameterMaps,
                       unsigned seed)
{
  for (auto& parameterMap: parameterMaps) {
    if (seed == 0) { // Use time
      parameterMap["RANDOM_SEED"].values = { (double) time(NULL) };
    } else {
      parameterMap["RANDOM_SEED"].values = { (double) seed};
    }

    unsigned numThreads = parameterMap.at("NUM_THREADS").dbl();
    unsigned numSimulations = parameterMap.at("NUM_SIMULATIONS").dbl();
    if (numThreads == 0) numThreads = numSimulations;
    unsigned simulationsRun = 0;

    // Construct symbol table
    SymbolTable symbolTable(numSimulations);
    for (auto& p: parameterMap.varyingParameters) {
      auto parameter = parameterMap.at(p);
      symbolTable.addSymbol(p, parameter.parentParameter, parameter.values);
    }

    // Needed for parameter printing
    unsigned parmPrint = (unsigned) parameterMap.at("PRINT_PARAMETERS").dbl();

    while(simulationsRun < numSimulations) {
      if (simulationsRun + numThreads > numSimulations)
        numThreads = numSimulations - simulationsRun;
      std::vector<std::thread> t(numThreads);

      for (unsigned i = 0; i < numThreads; ++i) {
        for (auto& symbol : symbolTable.table[simulationsRun + i]) {
          parameterMap[symbol.first].values = {symbol.second};
        }
        t[i] = std::thread(callSimulation, parameterMap,
                           parmPrint, simulationsRun + i);
      }

      for (unsigned i = 0; i < numThreads; ++i)  t[i].join();

      simulationsRun += numThreads;
    }
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

  TESTEQ(simulation.agents.size(), parameterMap.at("NUM_AGENTS").dbl(),
         successes, failures);

  simulation.agents[0]->age = 20;
  simulation.agents[0]->sex = MALE;
  simulation.agents[0]->desired_age = 50;
  simulation.agents[0]->sexual_orientation = HETEROSEXUAL;

  simulation.agents[1]->age = 30;
  simulation.agents[1]->sex = MALE;
  simulation.agents[1]->desired_age = 25;
  simulation.agents[1]->sexual_orientation = HOMOSEXUAL;

  simulation.agents[2]->age = 40;
  simulation.agents[2]->sex = FEMALE;
  simulation.agents[2]->desired_age = 30;
  simulation.agents[2]->sexual_orientation = HETEROSEXUAL;

  simulation.agents[3]->age = 25;
  simulation.agents[3]->sex = FEMALE;
  simulation.agents[3]->desired_age = 25;
  // simulation.agents[3]->relationship_length_factor = 0;
  simulation.agents[3]->sexual_orientation = HOMOSEXUAL;

  simulation.agents[4]->age = 30;
  simulation.agents[4]->sex = FEMALE;
  simulation.agents[4]->desired_age = 20;
  simulation.agents[4]->sexual_orientation = HOMOSEXUAL;

  simulation.agents[5]->age = 30;
  simulation.agents[5]->sex = MALE;
  simulation.agents[5]->desired_age = 20;
  simulation.agents[5]->sexual_orientation = HOMOSEXUAL;

  simulation.agents[6]->age = 25;
  simulation.agents[6]->sex = MALE;
  simulation.agents[6]->desired_age = 25;
  // simulation.agents[6]->relationship_length_factor = 0;
  simulation.agents[6]->sexual_orientation = HETEROSEXUAL;


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

  simulation.shapeRelationshipPeriod =
    CSVParser(parameterMap.at("SHAPE_REL_CSV").strValue.
              c_str(), ",", true).toDoubles();
  simulation.scaleRelationshipPeriod =
    CSVParser(parameterMap.at("SCALE_REL_CSV").strValue.
              c_str(), ",", true).toDoubles();
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

  parameterMap["END_DATE"].values = {2018};
  parameterMap["TIME_STEP"].values = {0.5};
  parameterMap["PRINT_PARAMETERS"].values = {0};
  parameterMap["ANALYZE_AFTER_INIT"].values = {0};
  parameterMap["ANALYZE_AT_END"].values = {0};
  parameterMap["OUTPUT_NUM_BREAKUPS"].values = {0};
  parameterMap["OUTPUT_NUM_MATINGPOOL"].values = {0};
  parameterMap["OUTPUT_TIMING_AFTER"].values = {0};
  Simulation s(parameterMap, 1);
  s.initializeAgents();
  s.calculateDemographics();
  TESTEQ(s.numInfectedMales > 0, true, successes, failures);
  TESTEQ(s.numInfectedFemales > 0, true, successes, failures);
  unsigned infected = s.numInfectedMales + s.numInfectedFemales;
  TESTEQ(infected < parameterMap.at("NUM_AGENTS").dbl(), true,
         successes, failures);
  s.simulate(false);
  s.calculateDemographics();
  TESTEQ(s.numInfectedMales + s.numInfectedFemales > infected, true,
         successes, failures);
  infected = s.numInfectedMales + s.numInfectedFemales;
  TESTEQ(infected < parameterMap.at("NUM_AGENTS").dbl(), true,
         successes, failures);

  printf("Successes: %u. Failures: %u.\n", successes, failures);
}
