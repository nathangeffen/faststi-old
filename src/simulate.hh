#ifndef SIMULATE_HH
#define SIMULATE_HH
#include <algorithm>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <cfloat>
#include <climits>
#include <cstdint>
#include <sys/time.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "CSV_Parser.hh"
#include "common.hh"
#include "linear.hh"
#include "agent.hh"
#include "parameters.hh"
#include "sample.hh"
#include "symboltable.hh"

double calcNumberSingles(DblMatrix, unsigned);
void setInitialInfection(Agent&,
                         const std::vector<double>&,
                         const std::vector<double>&,
                         const std::vector<double>&,
                         const std::vector<double>&);
void execSimulationSet(std::vector<ParameterMap>, unsigned);
void runTests(ParameterMap&);


class Simulation {
public:
  Simulation(const ParameterMap& parameter_map,
             const unsigned simulation_num) :
    parameterMap(parameter_map), simulationNum(simulation_num)
  {
    // Common to all simulations
    simulationName = parameterMap.at("SIMULATION_NAME").str();
    startDate = parameterMap.at("START_DATE").dbl();
    endDate = parameterMap.at("END_DATE").dbl();
    timeStep = parameterMap.at("TIME_STEP").dbl();
    currentDate = startDate;
    ageInterval = parameterMap.at("ANALYSIS_AGE_INTERVAL").dbl();
    failureThresholdScore = parameterMap.at("MATCH_SCORE_FAIL").dbl();
    poorThresholdScore = parameterMap.at("MATCH_SCORE_POOR").dbl();
    printNumBreakups = (parameterMap.at("OUTPUT_NUM_BREAKUPS").dbl() > 0.0)
      ? true : false;
    printNumMatings = (parameterMap.at("OUTPUT_NUM_MATINGPOOL").dbl() > 0.0)
      ? true : false;

    // Initialize thread local random number generator
    unsigned seed = parameterMap.at("RANDOM_SEED").dbl() * (simulationNum + 1);
    if (seed) rng.seed(seed);

    // Setup the events for this simulation
    setEvents();
    // Simulation class specific/
    initSimulation();
  }

  ~Simulation()
  {
    for (auto& agent: agents) delete agent;
  }

  virtual void initSimulation()
  {
    // Structures to track demographics and infections by age structure
    unsigned numIntervals = (MAX_AGE + 1) / ageInterval + 1;
    malesByAge = std::vector<unsigned>(numIntervals, 0);
    femalesByAge = std::vector<unsigned>(numIntervals, 0);
    infectedMalesByAge = std::vector<unsigned>(MAX_AGE + 1, 0);
    infectedFemalesByAge = std::vector<unsigned>(MAX_AGE + 1, 0);

    // Matching algorithm parameters
    neighbors = parameterMap.at("MATCH_NEIGHBORS").dbl();
    clusters = parameterMap.at("MATCH_CLUSTERS").dbl();

    // Infectiousness parameters
    het_male_infectiousness = parameterMap.at("HET_MALE_INFECTIOUSNESS").dbl();
    het_female_infectiousness =
      parameterMap.at("HET_FEMALE_INFECTIOUSNESS").dbl();
    hom_male_infectiousness = parameterMap.at("HOM_MALE_INFECTIOUSNESS").dbl();
    hom_female_infectiousness =
      parameterMap.at("HOM_FEMALE_INFECTIOUSNESS").dbl();

    // Single period parameters
    shapeSinglePeriodInitialMale =
      parameterMap.at("SHAPE_SINGLE_PERIOD_INITIAL_MALE").dbl();
    scaleSinglePeriodInitialMale =
      parameterMap.at("SCALE_SINGLE_PERIOD_INITIAL_MALE").dbl();
    shapeSinglePeriodInitialFemale =
      parameterMap.at("SHAPE_SINGLE_PERIOD_INITIAL_FEMALE").dbl();
    scaleSinglePeriodInitialFemale =
      parameterMap.at("SCALE_SINGLE_PERIOD_INITIAL_FEMALE").dbl();
    scaleSinglePeriodZeroDaysInitial =
      parameterMap.at("SCALE_SINGLE_PERIOD_ZERO_DAYS_INITIAL").dbl();
    scaleSinglePeriodZeroDaysDuring =
      parameterMap.at("SCALE_SINGLE_PERIOD_ZERO_DAYS_DURING").dbl();
    shapeSinglePeriodDuring =
      parameterMap.at("SHAPE_SINGLE_PERIOD_DURING").dbl();
    scaleSinglePeriodDuring =
      parameterMap.at("SCALE_SINGLE_PERIOD_DURING").dbl();
    meanSinglePeriodDeviation =
      parameterMap.at("MEAN_SINGLE_PERIOD").dbl();
    sdSinglePeriodDeviation =
      parameterMap.at("SD_SINGLE_PERIOD").dbl();

    // Relationship period parameters
    meanRelationshipPeriodDeviation =
      parameterMap.at("MEAN_RELATIONSHIP_PERIOD").dbl();
    sdRelationshipPeriodDeviation =
      parameterMap.at("SD_RELATIONSHIP_PERIOD").dbl();

    // Probability an agent is infected at initialisation if partner infected
    probInfectedIfPartnerInfected =
      parameterMap.at("PROB_INFECTED_IF_PARTNER").dbl();

    string s = parameterMap.at("DISTANCE_METHOD").str();
    if (s == "HEURISTIC")
      distanceMethod = HEURISTIC_DISTANCE;
    else if (s == "TABLE")
      distanceMethod = TABLE_DISTANCE;
    else {
      throw std::runtime_error ("Unknown distance method");
    }

    size_t max_p = (size_t) parameterMap.at("PARTNERS_RESERVE").dbl();
    if (max_p > 0) partnerships.partnerships.reserve(max_p);
    auto lf = parameterMap.at("PARTNERS_LOAD_FACTOR").dbl();
    if (lf > 0) partnerships.partnerships.max_load_factor(lf);

    // Age distribution matrices
    msmAgeDist = matrixFromCSV("MSM_AGE_DIST_CSV", ",", true);
    wswAgeDist = matrixFromCSV("WSW_AGE_DIST_CSV", ",", true);
    mswAgeDist = matrixFromCSV("MSW_AGE_DIST_CSV", ",", true);
    wsmAgeDist = matrixFromCSV("WSM_AGE_DIST_CSV", ",", true);

    // Probability single period is zero
    probZeroSinglePeriod = matrixFromCSV("PROB_ZERO_DAYS_SINGLE_CSV", ",", true);

    // Weibull parameters per age for setting relationship length
    shapeRelationshipPeriod = matrixFromCSV("SHAPE_REL_CSV", ",", true);
    scaleRelationshipPeriod = matrixFromCSV("SCALE_REL_CSV", ",", true);
  }

  void simulate(bool initAgents = true)
  {
    unsigned timing = parameterMap.at("OUTPUT_TIMING_DURING").dbl();

    struct timeval timeBegin, timeEnd;
    double elapsedTime;

    gettimeofday(&timeBegin, NULL);
    if (initAgents) initializeAgents();
    gettimeofday(&timeEnd, NULL);
    elapsedTime = timeEnd.tv_sec - timeBegin.tv_sec;
    if (parameterMap.at("OUTPUT_INIT_TIMING").dbl()) {
      printf("%s,TIMING,INIT,%u,%f\n", simulationName.c_str(),
             simulationNum, elapsedTime);
    }
    unsigned outputAgents = parameterMap.at("OUTPUT_AGENTS_AFTER_INIT").dbl();
    if (outputAgents) printAgents(agents, simulationNum, startDate, stdout);
    /********************/
    /* Count WSW */
    unsigned wswcount = 0;
    for (auto& a: agents) {
      if (a->sex == FEMALE && a->sexual_orientation == HOMOSEXUAL) ++wswcount;
    }
    /*******************/
    unsigned analyzeAgents = parameterMap.at("ANALYZE_AFTER_INIT").dbl();
    if (analyzeAgents) analysis();
    /* Main loop */

    // Make sure main loop uses integer arithmetic, rather than floats
    // though probably makes little diff.
    unsigned num_iterations = (endDate - startDate) / timeStep;
    outputAgents = parameterMap.at("OUTPUT_AGENTS_DURING_SIM").dbl(0);
    unsigned outputFrequency =
      parameterMap.at("OUTPUT_AGENTS_DURING_SIM").dbl(1);
    unsigned analyzeFrequency = parameterMap.at("ANALYZE_DURING_SIM").dbl();
    for (unsigned i = 0; i < num_iterations; ++i, currentDate += timeStep) {
      for (auto& e: events) e(this);
      if (timing > 0 && (i + 1) % timing == 0) {
        gettimeofday(&timeEnd, NULL);
        elapsedTime = timeEnd.tv_sec - timeBegin.tv_sec;
        printf("%s,TIMING,%u,%u,%f\n",
               simulationName.c_str(), i, simulationNum, elapsedTime);
      }
      if (outputAgents > 0 && (i + 1) % outputFrequency == 0)
        printAgents(agents, simulationNum, currentDate, stdout);
      if (analyzeFrequency && (i + 1) % analyzeFrequency == 0) analysis();
    }

    gettimeofday(&timeEnd, NULL);
    elapsedTime = timeEnd.tv_sec - timeBegin.tv_sec;
    if (parameterMap.at("OUTPUT_TIMING_AFTER").dbl()) {
      printf("%s,TIMING,AFTER,%u,%f\n",
             simulationName.c_str(), simulationNum, elapsedTime);
    }

    /* Wrap up */
    outputAgents = parameterMap.at("OUTPUT_AGENTS_AT_END").dbl();
    if (outputAgents) printAgents(agents, simulationNum, endDate, stdout);
    analyzeAgents = parameterMap.at("ANALYZE_AT_END").dbl();
    if (analyzeAgents) analysis();
  }


  inline DblMatrix matrixFromCSV(const char* key,
                                 const char* delim,
                                 const bool hasHeader)
  {
    std::string filename = parameterMap.at(key).str();
    CSVParser csvParser(filename.c_str(), delim, hasHeader);
    DblMatrix result = csvParser.toDoubles();
    return result;
  }

  virtual void initializeAgents()
  {
    // Read CSV files of agent demographics and sexual behaviour
    // into matrices of doubles
    DblMatrix demographics = matrixFromCSV("AGENT_DATA_CSV", ";", true);
    assert(demographics.size());
    DblMatrix singles = matrixFromCSV("SINGLES_DATA_CSV", ",", true);
    DblMatrix partners = matrixFromCSV("PARTNERS_DATA_CSV", ",", true);
    DblMatrix mm = matrixFromCSV("MSM_DATA_CSV", ";", false);
    DblMatrix ww = matrixFromCSV("WSW_DATA_CSV", ";", false);
    DblMatrix mw = matrixFromCSV("MSW_DATA_CSV", ";", false);
    DblMatrix wm = matrixFromCSV("WSM_DATA_CSV", ";", false);

    assert(singles.size());

    // Calculate initial infection rates matrices
    DblMatrix initialRates = matrixFromCSV("INITIAL_INFECTION_RATES_CSV",
                                           ",", true);
    double scaleInitialRates = parameterMap.at("SCALE_INITIAL_RATES").dbl();
    std::vector<double> initialInfectionRatesMSW(MAX_AGE + 1, 0.0);
    std::vector<double> initialInfectionRatesMSM(MAX_AGE + 1, 0.0);
    std::vector<double> initialInfectionRatesWSM(MAX_AGE + 1, 0.0);
    std::vector<double> initialInfectionRatesWSW(MAX_AGE + 1, 0.0);
    {
      unsigned i = 0;
      for (auto& row : initialRates) {
        for (;i <= row[0] && i < MAX_AGE + 1; ++i) {
          assert(i <= MAX_AGE);
          initialInfectionRatesMSW[i] = row[1] * scaleInitialRates;
          initialInfectionRatesMSM[i] = row[2] * scaleInitialRates;
          initialInfectionRatesWSM[i] = row[3] * scaleInitialRates;
          initialInfectionRatesWSW[i] = row[4] * scaleInitialRates;
        }
      }
      for (; i <= MAX_AGE; ++i) {
        initialInfectionRatesMSW[i] = 0.0;
        initialInfectionRatesMSM[i] = 0.0;
        initialInfectionRatesWSM[i] = 0.0;
        initialInfectionRatesWSW[i] = 0.0;
      }
    }

    // Calculate number of agents in relationships and number that are single
    unsigned X = parameterMap.at("NUM_AGENTS").values[0];
    unsigned S = calcNumberSingles(demographics, X);
    agents.reserve(X);

    // Get age structure, sex & sexual orientation information for single agents
    std::vector<double> ageRange;
    for (unsigned i = 12; i <= 100; ++i) ageRange.push_back(i);
    auto ageShare = getCol(singles, 1);
    auto femRatio = getCol(singles, 2);
    auto msmRate = getCol(singles, 3);
    auto wswRate = getCol(singles, 4);

    // Relationship period parameter for beginning of period
    scaleModifierRelationshipPeriod =
      parameterMap.at("SCALE_RELATIONSHIP_PERIOD_INITIAL").dbl();

    // Create the single agents
    createAgents(agents,
                 0, S,
                 ageRange, ageShare, femRatio,
                 wswRate, msmRate, ww, mw, wm, mm,
                 initialInfectionRatesMSW, initialInfectionRatesMSM,
                 initialInfectionRatesWSM, initialInfectionRatesWSW, false);

    // Get age structure, sex & sexual orientation information for paired agents
    ageShare = femRatio = msmRate = wswRate = {};
    ageShare = getCol(partners, 1);
    femRatio = getCol(partners, 2);
    msmRate = getCol(partners, 3);
    wswRate = getCol(partners, 4);

    // Create the paired agents and their partners
    createAgents(agents,
                 S, X,
                 ageRange, ageShare, femRatio,
                 wswRate, msmRate, ww, mw, wm, mm,
                 initialInfectionRatesMSW, initialInfectionRatesMSM,
                 initialInfectionRatesWSM, initialInfectionRatesWSW, true);

    // Optimise space
    agents.shrink_to_fit();

    // Get the initial demographic information for reporting
    calculateDemographics();

    // Relationship period parameter for during simulation
    scaleModifierRelationshipPeriod =
      parameterMap.at("SCALE_RELATIONSHIP_PERIOD_DURING").dbl();
  }

  inline void initAgent(Agent *agent,
                        const unsigned id,
                        Sample& sampleAgeshare,
                        const std::vector<double>& femRatio,
                        const std::vector<double>& wswRate,
                        const std::vector<double>& msmRate,
                        const std::vector<double>& initialInfectionRatesMSW,
                        const std::vector<double>& initialInfectionRatesMSM,
                        const std::vector<double>& initialInfectionRatesWSM,
                        const std::vector<double>& initialInfectionRatesWSW,
                        std::vector<Sample>& sample_matWW,
                        std::vector<Sample>& sample_matMW,
                        std::vector<Sample>& sample_matWM,
                        std::vector<Sample>& sample_matMM)
  {
    std::uniform_real_distribution<double> uni;

    agent->id = id;

    // Age
    unsigned age = sampleAgeshare() + 12;
    agent->age = age;
    // Sex
    unsigned sex = uni(rng) < femRatio[age - 12] ? FEMALE : MALE;
    agent->sex = sex;
    // Sexual_Orientation
    unsigned sexual_orientation;
    if (sex == FEMALE) {
      sexual_orientation = uni(rng) < wswRate[age - 12]
                                      ? HOMOSEXUAL : HETEROSEXUAL;
    } else {
      sexual_orientation = uni(rng) < msmRate[age - 12]
                                      ? HOMOSEXUAL : HETEROSEXUAL;
    }
    agent->sexual_orientation = sexual_orientation;
    setInitialInfection(*agent, initialInfectionRatesMSW,
                        initialInfectionRatesMSM, initialInfectionRatesWSM,
                        initialInfectionRatesWSW);
    // Desired age of partner
    if (sex == FEMALE && sexual_orientation == HOMOSEXUAL)
      agent->desired_age  = sample_matWW[age - 12]() + 12;
    else if (sex == FEMALE && sexual_orientation == HETEROSEXUAL)
      agent->desired_age = sample_matWM[age - 12]() + 12;
    else if (sex == MALE && sexual_orientation == HETEROSEXUAL)
      agent->desired_age = sample_matMW[age - 12]() + 12;
    else
      agent->desired_age = sample_matMM[age - 12]() + 12;
    // Agent's deviation from mean for period spent in relationships or single
    std::normal_distribution<double>
      normSingle(meanSinglePeriodDeviation, sdSinglePeriodDeviation);
    std::normal_distribution<double>
      normRelationship(meanRelationshipPeriodDeviation,
                       sdRelationshipPeriodDeviation);
    agent->singlePeriodDeviation = normSingle(rng);
    agent->relationshipPeriodDeviation = normRelationship(rng);
  }

  inline void conditionalInfectPartner(Agent *agent)
  {
    std::uniform_real_distribution<double> uni;
    if (uni(rng) < probInfectedIfPartnerInfected) agent->infected = true;
  }

  void createAgents(AgentVector& agents,
                    const unsigned fromAgent,
                    const unsigned toAgent,
                    const std::vector<double>& ageRange,
                    const std::vector<double>& ageShare,
                    const std::vector<double>& femRatio,
                    const std::vector<double>& wswRate,
                    const std::vector<double>& msmRate,
                    const DblMatrix& matWW,
                    const DblMatrix& matMW,
                    const DblMatrix& matWM,
                    const DblMatrix& matMM,
                    const std::vector<double>& initialInfectionRatesMSW,
                    const std::vector<double>& initialInfectionRatesMSM,
                    const std::vector<double>& initialInfectionRatesWSM,
                    const std::vector<double>& initialInfectionRatesWSW,
                    bool initial_relation)
  {
    std::uniform_real_distribution<double> uni;
    Sample sample_ageshare(ageShare, &rng);
    vector<Sample> sample_matWW(matWW[0].size());
    vector<Sample> sample_matMW(matMW[0].size());
    vector<Sample> sample_matWM(matWM[0].size());
    vector<Sample> sample_matMM(matMM[0].size());
    std::vector<double> placeholder(matMM.size(), 0.000001);

    for (unsigned i = 0; i < matWW[0].size(); ++i) {
      sample_matWW[i].init(getCol(matWW,i), &rng);
      sample_matWM[i].init(getCol(matWM,i), &rng);
      sample_matMW[i].init(getCol(matMW,i), &rng);
      sample_matMM[i].init(getCol(matMM,i), &rng);
    }


    unsigned increment = initial_relation ? 2 : 1;
    for (unsigned i = fromAgent; (i + increment - 1) < toAgent; i += increment) {
      Agent *agent = new Agent();
      initAgent(agent, i + 1,
                sample_ageshare,
                femRatio,
                wswRate,
                msmRate,
                initialInfectionRatesMSW,
                initialInfectionRatesMSM,
                initialInfectionRatesWSM,
                initialInfectionRatesWSW,
                sample_matWW,
                sample_matMW,
                sample_matWM,
                sample_matMM);

      if (initial_relation) {
        // Relationship
        agent->initial_relationship = true;
        Agent* partner = new Agent();
        initAgent(partner, i + 2,
                  sample_ageshare,
                  femRatio,
                  wswRate,
                  msmRate,
                  initialInfectionRatesMSW,
                  initialInfectionRatesMSM,
                  initialInfectionRatesWSM,
                  initialInfectionRatesWSW,
                  sample_matWW,
                  sample_matMW,
                  sample_matWM,
                  sample_matMM);
        // Orientation = partner's orientation
        partner->sexual_orientation = agent->sexual_orientation;
        // Partner sex
        if (agent->sexual_orientation == HETEROSEXUAL) {
          if (agent->sex == MALE)
            partner->sex = FEMALE;
          else
            partner->sex = MALE;
        } else {
          partner->sex = agent->sex;
        }
        partner->age = agent->desired_age;
        // Partner in relationship
        partner->initial_relationship = true;
        // Preferred age of partner
        partner->desired_age = agent->age;
        // partner infection risk parameters
        setInitialInfection(*partner, initialInfectionRatesMSW,
                            initialInfectionRatesMSM, initialInfectionRatesWSM,
                            initialInfectionRatesWSW);

        makePartner(agent, partner, distance(agent, partner));
        // Correct relationship time because this is in the middle of relationship
        std::uniform_real_distribution<double>
          uni2(currentDate, agent->relationshipChangeDate);
        agent->relationshipChangeDate = uni2(rng);
        partner->relationshipChangeDate = agent->relationshipChangeDate;
        if (agent->infected == true && partner->infected == false) {
          conditionalInfectPartner(partner);
        } else if (agent->infected == false && partner->infected == true) {
          conditionalInfectPartner(agent);
        }
      } else {
        if (agent->sex == MALE) {
          agent->setSinglePeriod(currentDate,
                                 shapeSinglePeriodInitialMale,
                                 scaleSinglePeriodInitialMale,
                                 probZeroSinglePeriod,
                                 scaleSinglePeriodZeroDaysInitial);
        } else {
          agent->setSinglePeriod(currentDate,
                                 shapeSinglePeriodInitialFemale,
                                 scaleSinglePeriodInitialFemale,
                                 probZeroSinglePeriod,
                                 scaleSinglePeriodZeroDaysInitial);
        }
      }
      agents.push_back(agent);
      if (agent->partner) agents.push_back(agent->partner);
    }
  }

  void calculateDemographics()
  {
    for (unsigned i = 0; i < agents.size(); ++i) {
      unsigned age = agents[i]->age;
      if (agents[i]->sex == MALE) {
        ++numMales;
        ++malesByAge[age / ageInterval];
        if (agents[i]->sexual_orientation == HETEROSEXUAL)
          ++numMsw;
        else
          ++numMsm;
        if (agents[i]->infected) {
          ++numInfectedMales;
          ++infectedMalesByAge[age / ageInterval];
          if (agents[i]->sexual_orientation == HETEROSEXUAL)
            ++numInfectedMsw;
          else
            ++numInfectedMsm;
        }
      } else {
        ++numFemales;
        ++femalesByAge[age /  ageInterval];
        if (agents[i]->sexual_orientation == HETEROSEXUAL)
          ++numWsm;
        else
          ++numWsw;
        if (agents[i]->infected) {
          ++numInfectedFemales;
          ++infectedFemalesByAge[age / ageInterval];
          if (agents[i]->sexual_orientation == HETEROSEXUAL)
            ++numInfectedWsm;
          else
            ++numInfectedWsw;
        }
      }
    }
  }

  void trackRiskFactors(Agent* agent)
  {
    if (agent->sex == MALE) {
      ++numInfectedMales;
      ++infectedMalesByAge[ (int) agent->age / ageInterval];
      if (agent->sexual_orientation == HETEROSEXUAL)
        ++numInfectedMsw;
      else
        ++numInfectedMsm;
    } else {
      ++numInfectedFemales;
      ++infectedFemalesByAge[ (int) agent->age / ageInterval];
      if (agent->sexual_orientation == HETEROSEXUAL)
        ++numInfectedWsm;
      else
        ++numInfectedWsw;
    }
  }

  AgentVector getUnmatchedAgents()
  {
    AgentVector matingPool;
    for (auto& agent : agents) {
      if (agent->isMatchable(currentDate))
        matingPool.push_back(agent);
    }
    return matingPool;
  }

  AgentVector getShuffledUnmatchedAgents()
  {
    AgentVector matingPool = getUnmatchedAgents();
    shuffle(matingPool.begin(), matingPool.end(), rng);
    // Remove back agent if odd
    if (matingPool.size() % 2 == 1) matingPool.pop_back();
    if (printNumMatings) {
      printf("%s,MATINGPOOL,,%u,%.3f,%lu\n", simulationName.c_str(),
             simulationNum, currentDate, matingPool.size());
    }
    return matingPool;
  }

  PartnershipScore
  closestPairMatchN(std::vector<Agent *>::iterator from,
		    std::vector<Agent *>::iterator to)
  {
    double smallest_val = DBL_MAX;
    std::vector<Agent *>::iterator closest_agent = to;

    for (auto it = from + 1; it != to; ++it) {
      if ( (*it)->partner == NULL) {
        double d = distance(*from, *it);
        if (d < smallest_val) {
          smallest_val = d;
          closest_agent = it;
        }
      }
    }
    PartnershipScore partnershipScore = {closest_agent, smallest_val};
    return partnershipScore;
  }


  void analysis()
  {
    double prevalence = (double) (numInfectedMales + numInfectedFemales) /
      (numMales + numFemales);
    double malePrevalence = (double) numInfectedMales / numMales;
    double femalePrevalence = (double) numInfectedFemales / numFemales;
    double msmPrevalence = (double) numInfectedMsm / numMsm;
    double wswPrevalence = (double) numInfectedWsw / numWsw;

    printf("%s,ANALYSIS,INFECTED,%d,%.3f,%u\n",
           simulationName.c_str(), simulationNum, currentDate,
           numInfectedMales + numInfectedFemales);
    printf("%s,ANALYSIS,PREVALENCE,%d,%.3f,%f\n",
           simulationName.c_str(), simulationNum, currentDate, prevalence);
    printf("%s,ANALYSIS,MALEPREVALENCE,%d,%.3f,%f\n",
           simulationName.c_str(),simulationNum, currentDate, malePrevalence);
    printf("%s,ANALYSIS,FEMALEPREVALENCE,%d,%.3f,%f\n",
           simulationName.c_str(),simulationNum, currentDate, femalePrevalence);
    printf("%s,ANALYSIS,MSMPREVALENCE,%d,%.3f,%f\n",
           simulationName.c_str(),simulationNum, currentDate, msmPrevalence);
    printf("%s,ANALYSIS,WSWPREVALENCE,%d,%.3f,%f\n",
           simulationName.c_str(),simulationNum, currentDate, wswPrevalence);
    printf("%s,ANALYSIS,PARTNERSHIPS,%d,%.3f,%u\n",
           simulationName.c_str(),simulationNum, currentDate, totalPartnerships);
    printf("%s,ANALYSIS,MSMPARTNERSHIPS,%d,%.3f,%u\n",
           simulationName.c_str(),simulationNum, currentDate,
           totalMsmPartnerships);
    printf("%s,ANALYSIS,WSWPARTNERSHIPS,%d,%.3f,%u\n",
           simulationName.c_str(), simulationNum, currentDate,
           totalWswPartnerships);

    for (unsigned i = 0; i < malesByAge.size(); ++i) {
      ostringstream ssmale, ssfemale;
      if (malesByAge[i] == 0)
        ssmale << "NA";
      else
        ssmale << std::setprecision(3)
               << ( (double) infectedMalesByAge[i] / malesByAge[i] );
      if (femalesByAge[i] == 0)
        ssfemale << "NA";
      else
        ssfemale << std::setprecision(3)
                 << ((double) infectedFemalesByAge[i] / femalesByAge[i]);
      printf("%s,ANALYSIS,MALE_AGE_%03u-%03u,%d,%.3f,%s\n",
             simulationName.c_str(),
             i * ageInterval, i * ageInterval + ageInterval - 1,
             simulationNum, currentDate, ssmale.str().c_str());
      printf("%s,ANALYSIS,FEMALE_AGE_%03u-%03u,%d,%.3f,%s\n",
             simulationName.c_str(),
             i * ageInterval, i * ageInterval + ageInterval - 1,
             simulationNum, currentDate, ssfemale.str().c_str());
    }
    printf("%s,ANALYSIS,SCORE,%d,%.3f,%f\n",
           simulationName.c_str(), simulationNum, currentDate,
           totalPartnershipScore / totalPartnerships);
    printf("%s,ANALYSIS,FAILED,%d,%.3f,%u\n",
           simulationName.c_str(), simulationNum, currentDate, failedMatches);
    printf("%s,ANALYSIS,POOR,%d,%.3f,%u\n",
           simulationName.c_str(), simulationNum, currentDate, poorMatches);
  }


  void makePartner(Agent* a, Agent *b,
                   const double score)
  {
    assert(a->partner == NULL);
    assert(b->partner == NULL);

    if (score < failureThresholdScore) {
      if (score > poorThresholdScore) ++poorMatches;
      totalPartnershipScore += score;
      partnerships.insert(a->id, b->id);
      a->partner = b;
      b->partner = a;
      a->setRelationshipPeriod(currentDate, shapeRelationshipPeriod,
                               scaleRelationshipPeriod,
                               scaleModifierRelationshipPeriod);
      ++a->num_partners;
      ++b->num_partners;
      ++totalPartnerships;
      if (a->sex == b->sex) {
        if (a->sex == MALE)
          ++totalMsmPartnerships;
        else
          ++totalWswPartnerships;
      } else {
        ++totalMswPartnerships;
      }
    } else {
      ++failedMatches;
    }
  }


  double distance(const Agent *a, const Agent *b) const
  {
    return distanceMethod == HEURISTIC_DISTANCE
      ? heuristicDistance(a, b) : tableDistance(a, b);
  }

  double heuristicDistance(const Agent *a, const Agent *b) const
  {
    double score = 0.0;

    score +=  (fabs(a->desired_age - b->age) +
               fabs(b->desired_age - a->age)) / 2.0;
    if (a->sexual_orientation != b->sexual_orientation) {
      score += 50.0;
    } else if (a->sexual_orientation == HETEROSEXUAL) {
      if (a->sex == b->sex)
        score += 50.0;
    } else if (a->sex != b->sex) {
      score += 50.0;
    }
    if (partnerships.exists(a->id, b->id))
      score += 50.0;
    return score;
  }

  double tableDistance(const Agent *a, const Agent *b) const
  {
    double score = 0.0;
    unsigned a_age = std::min( (unsigned) a->age, (unsigned) MAX_AGE) - MIN_AGE;
    unsigned b_age = std::min( (unsigned) b->age,  (unsigned) MAX_AGE) - MIN_AGE;

    if (a->sex == MALE and b->sex == FEMALE) {
      score += (mswAgeDist[a_age][b_age] + wsmAgeDist[b_age][a_age]) * 25;
    } else if (a->sex == FEMALE and b->sex == MALE) {
      score += (wsmAgeDist[a_age][b_age] + mswAgeDist[b_age][a_age]) * 25;
    } else if (a->sex == MALE and b->sex == MALE) {
      score += (msmAgeDist[a_age][b_age] + msmAgeDist[b_age][a_age]) * 25;
    } else if (a->sex == FEMALE and b->sex == FEMALE) {
      score += (wswAgeDist[a_age][b_age] + wswAgeDist[b_age][a_age]) * 25;
    } else {
      throw std::runtime_error("Unknown sex combination");
    }

    if (a->sexual_orientation != b->sexual_orientation) {
      score += 50.0;
    } else if (a->sexual_orientation == HETEROSEXUAL) {
      if (a->sex == b->sex)
        score += 50.0;
    } else if (a->sex != b->sex) {
      score += 50.0;
    }

    if (partnerships.exists(a->id, b->id))
      score += 50.0;

    return score;
  }

  double clusterValue(const Agent *a) const
  {
    return ( (double) a->desired_age / 100.0 + a->age / 100.0) / 2.0
      + a->sexual_orientation;
  }

  virtual void setEvents();

  std::string simulationName;
  AgentVector agents;
  Partnerships partnerships;
  ParameterMap parameterMap;
  unsigned simulationNum;
  unsigned ageInterval;
  unsigned distanceMethod;
  double startDate;
  double endDate;
  double timeStep;

  double shapeSinglePeriodInitialMale;
  double scaleSinglePeriodInitialMale;
  double shapeSinglePeriodInitialFemale;
  double scaleSinglePeriodInitialFemale;
  double shapeSinglePeriodDuring;
  double scaleSinglePeriodDuring;
  double meanSinglePeriodDeviation;
  double sdSinglePeriodDeviation;

  double scaleModifierRelationshipPeriod;
  double meanRelationshipPeriodDeviation;
  double sdRelationshipPeriodDeviation;

  double het_male_infectiousness;
  double hom_male_infectiousness;
  double het_female_infectiousness;
  double hom_female_infectiousness;
  double probInfectedIfPartnerInfected;

  bool printNumMatings;
  bool printNumBreakups;

  double currentDate;
  double failureThresholdScore;
  double poorThresholdScore;
  double totalPartnershipScore = 0.0;
  unsigned clusters;
  unsigned neighbors;
  unsigned totalBreakups = 0;
  unsigned totalPartnerships = 0;
  unsigned totalMswPartnerships = 0;
  unsigned totalMsmPartnerships = 0;
  unsigned totalWswPartnerships = 0;
  unsigned numMales = 0;
  unsigned numFemales = 0;
  unsigned numMsm = 0;
  unsigned numMsw = 0;
  unsigned numWsm = 0;
  unsigned numWsw = 0;
  unsigned numInfectedMales = 0;
  unsigned numInfectedFemales = 0;
  unsigned numInfectedMsm = 0;
  unsigned numInfectedMsw = 0;
  unsigned numInfectedWsm = 0;
  unsigned numInfectedWsw = 0;
  unsigned numShortBreakShortPartnership = 0;
  unsigned numShortBreakLongPartnership = 0;
  unsigned numLongBreakShortPartnership = 0;
  unsigned numLongBreakLongPartnership = 0;
  unsigned numInfectedShortBreakShortPartnership = 0;
  unsigned numInfectedShortBreakLongPartnership = 0;
  unsigned numInfectedLongBreakShortPartnership = 0;
  unsigned numInfectedLongBreakLongPartnership = 0;
  unsigned poorMatches = 0;
  unsigned failedMatches = 0;
  std::vector<unsigned> malesByAge;
  std::vector<unsigned> femalesByAge;
  std::vector<unsigned> infectedMalesByAge;
  std::vector<unsigned> infectedFemalesByAge;

  DblMatrix shapeRelationshipPeriod;
  DblMatrix scaleRelationshipPeriod;
  DblMatrix mswAgeDist;
  DblMatrix wsmAgeDist;
  DblMatrix msmAgeDist;
  DblMatrix wswAgeDist;
  DblMatrix probZeroSinglePeriod;
  double scaleSinglePeriodZeroDaysInitial;
  double scaleSinglePeriodZeroDaysDuring;

  std::vector< std::function<void(Simulation*)> > events;
};

#endif
