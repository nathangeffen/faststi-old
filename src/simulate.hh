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
#include <cstring>
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

double calcNumberSingles(const DblMatrix&, const unsigned);
void setInitialInfection(Agent&,
                         const std::vector<double>&,
                         const std::vector<double>&,
                         const std::vector<double>&,
                         const std::vector<double>&);
void execSimulationSet(std::vector<ParameterMap>, unsigned);
void runTests(ParameterMap&);

void writeCsvLine(const std::string&,
                  const unsigned,
                  const double,
                  const std::string&,
                  const std::string&,
                  const std::ostringstream& = std::ostringstream(""));

/**
   Used to track number of males and females, and number of infected of these,
   by age interval.
 */

struct AgeStatistics {
  AgeStatistics() {
    const size_t SIZE_TO_FREE = (unsigned) NUM_INTERVALS * sizeof(unsigned);
    memset(malesByAge, 0, SIZE_TO_FREE);
    memset(femalesByAge, 0, SIZE_TO_FREE);
    memset(infectedMalesByAge, 0, SIZE_TO_FREE);
    memset(infectedFemalesByAge, 0, SIZE_TO_FREE);
  }
  unsigned malesByAge[NUM_INTERVALS];
  unsigned femalesByAge[NUM_INTERVALS];
  unsigned infectedMalesByAge[NUM_INTERVALS];
  unsigned infectedFemalesByAge[NUM_INTERVALS];
};

/**
   Simulation engine for Faststi.
 */

class Simulation {

public:
  std::string simulationName;
  AgentVector agents;
  Partnerships partnerships;
  ParameterMap parameterMap;
  unsigned simulationNum;
  unsigned distanceMethod;
  double startDate;
  double endDate;
  double timeStep;

  double scaleSinglePeriodZeroDaysInitial;
  double scaleSinglePeriodZeroDaysDuring;
  double scaleSinglePeriodInitial;
  double scaleSinglePeriodDuring;

  double meanSinglePeriodDeviation;
  double sdSinglePeriodDeviation;

  double scaleModifierRelationshipPeriod;
  double meanRelationshipPeriodDeviation;
  double sdRelationshipPeriodDeviation;

  double hetMaleInfectiousness;
  double homMaleInfectiousness;
  double hetFemaleInfectiousness;
  double homFemaleInfectiousness;
  double probInfectedIfPartnerInfected;

  bool printNumMatings;
  bool printNumBreakups;
  bool analyzeTruncatedAgeInit = false;
  bool analyzeTruncatedAgeDuring = false;
  bool analyzeTruncatedAgeAfter = false;
  bool incTruncatedAge = false;

  double currentDate;
  double failureThresholdScore;
  double poorThresholdScore;
  double totalPartnershipScore = 0.0;

  unsigned stabilizationSteps;
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
  unsigned minAge;
  unsigned maxAge;

  DblMatrix weibullSinglePeriodInitial;
  DblMatrix weibullSinglePeriodDuring;
  DblMatrix shapeRelationshipPeriod;
  DblMatrix scaleRelationshipPeriod;
  DblMatrix mswAgeDist;
  DblMatrix wsmAgeDist;
  DblMatrix msmAgeDist;
  DblMatrix wswAgeDist;
  DblMatrix probZeroSinglePeriod;

  std::vector< std::function<void(Simulation*)> > events;

  /**
      Constructor for simulation engine.

      @param parameter_map[in] User defined parameters for this simulation.
      @param simulation_num[in] Unique simulation number.
   */
  explicit Simulation(const ParameterMap& parameter_map,
                      const unsigned simulation_num) :
    parameterMap(parameter_map), simulationNum(simulation_num)
  {
    // Common to all simulations
    simulationName = parameterMap.at("SIMULATION_NAME").str();
    startDate = parameterMap.at("START_DATE").dbl();
    endDate = parameterMap.at("END_DATE").dbl();
    timeStep = parameterMap.at("TIME_STEP").dbl();
    stabilizationSteps = parameterMap.at("STABILIZATION_STEPS").dbl();
    currentDate = startDate - (double) stabilizationSteps * timeStep;
    failureThresholdScore = parameterMap.at("MATCH_SCORE_FAIL").dbl();
    poorThresholdScore = parameterMap.at("MATCH_SCORE_POOR").dbl();
    printNumBreakups = (parameterMap.at("OUTPUT_NUM_BREAKUPS").dbl() > 0.0)
      ? true : false;
    printNumMatings = (parameterMap.at("OUTPUT_NUM_MATINGPOOL").dbl() > 0.0)
      ? true : false;
    analyzeTruncatedAgeInit = parameterMap.at("ANALYZE_TRUNCATED_AGE_INIT").isSet();
    analyzeTruncatedAgeDuring = parameterMap.at("ANALYZE_TRUNCATED_AGE_DURING").isSet();
    analyzeTruncatedAgeAfter = parameterMap.at("ANALYZE_TRUNCATED_AGE_AFTER").isSet();


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

  /**
     Initializes engine's class members with user defined parameters.
   */

  virtual void initSimulation()
  {
    // Structures to track demographics and infections by age structure

    // Matching algorithm parameters
    neighbors = parameterMap.at("MATCH_NEIGHBORS").dbl();
    clusters = parameterMap.at("MATCH_CLUSTERS").dbl();

    // Infectiousness parameters
    hetMaleInfectiousness = parameterMap.at("HET_MALE_INFECTIOUSNESS").dbl();
    hetFemaleInfectiousness =
      parameterMap.at("HET_FEMALE_INFECTIOUSNESS").dbl();
    homMaleInfectiousness = parameterMap.at("HOM_MALE_INFECTIOUSNESS").dbl();
    homFemaleInfectiousness =
      parameterMap.at("HOM_FEMALE_INFECTIOUSNESS").dbl();

    // Single period parameters
    scaleSinglePeriodZeroDaysInitial =
      parameterMap.at("SCALE_SINGLE_PERIOD_ZERO_DAYS_INITIAL").dbl();
    scaleSinglePeriodZeroDaysDuring =
      parameterMap.at("SCALE_SINGLE_PERIOD_ZERO_DAYS_DURING").dbl();
    scaleSinglePeriodInitial =
      parameterMap.at("SCALE_SINGLE_PERIOD_INITIAL").dbl();
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

    // Relationship period parameter for beginning of period
    scaleModifierRelationshipPeriod =
      parameterMap.at("SCALE_RELATIONSHIP_PERIOD_INITIAL").dbl();

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

    // Minimum and maximum ages
    minAge = parameterMap.at("MIN_AGE_TRUNCATE").dbl();
    if (minAge == 0) minAge = MIN_AGE;
    maxAge = parameterMap.at("MAX_AGE_TRUNCATE").dbl();
    if (maxAge == 0) maxAge = MAX_AGE;
    incTruncatedAge = parameterMap.at("INC_TRUNCATED_AGE").isSet();

    // Age distribution matrices
    msmAgeDist = matrixFromCSV("MSM_AGE_DIST_CSV", ",", true);
    wswAgeDist = matrixFromCSV("WSW_AGE_DIST_CSV", ",", true);
    mswAgeDist = matrixFromCSV("MSW_AGE_DIST_CSV", ",", true);
    wsmAgeDist = matrixFromCSV("WSM_AGE_DIST_CSV", ",", true);
    // size_t n = maxAge - MIN_AGE;
    // truncateMatrix(msmAgeDist, n, n);
    // truncateMatrix(wswAgeDist, n, n);
    // truncateMatrix(mswAgeDist, n, n);
    // truncateMatrix(wsmAgeDist, n, n);
    // Probability single period is zero
    probZeroSinglePeriod = matrixFromCSV("PROB_ZERO_DAYS_SINGLE_CSV", ",", true);

    // Weibull parameters per age for setting single period
    weibullSinglePeriodInitial =
      matrixFromCSV("WEIBULL_SINGLE_INITIAL_CSV", ",", true);

    weibullSinglePeriodDuring =
      matrixFromCSV("WEIBULL_SINGLE_DURING_CSV", ",", true);


    // Weibull parameters per age for setting relationship length
    shapeRelationshipPeriod = matrixFromCSV("SHAPE_REL_CSV", ",", true);
    scaleRelationshipPeriod = matrixFromCSV("SCALE_REL_CSV", ",", true);
  }

  /**
     Engine of the simulation. It iterates by the sizeof the time step from
     the beginning date of the simulation till the end date, executing events
     on each iteration.

     @param initAgents[in] Indicates whether or not to execute the agent
     initialization virtual function.
   */
  void simulate(const bool initAgents = true)
  {
    unsigned timing = parameterMap.at("OUTPUT_TIMING_DURING").dbl();

    struct timeval timeBegin, timeEnd;
    double elapsedTime;

    gettimeofday(&timeBegin, NULL);
    if (initAgents) initializeAgents();
    gettimeofday(&timeEnd, NULL);
    elapsedTime = timeEnd.tv_sec - timeBegin.tv_sec;
    if (parameterMap.at("OUTPUT_INIT_TIMING").dbl()) {
      csvout("TIMING", "INIT", elapsedTime);
    }
    if (parameterMap.at("OUTPUT_AGENTS_INIT").isSet()) {
      printAgents(agents, simulationNum, startDate);
    }
    if (parameterMap.at("ANALYZE_INIT").isSet()) {
      analysis(true, true, false, analyzeTruncatedAgeInit);
    }

    /* Main loop */
    // Make sure main loop uses integer arithmetic, rather than floats
    // though probably makes little diff.
    unsigned numIterations =
      ( (endDate - startDate) / timeStep ) + stabilizationSteps;
#ifdef DEBUG
    if (parameterMap.at("PRINT_NUM_ITERATIONS_AND_EXIT").isSet()) {
      std::cerr << "Number of iterations: " << numIterations << std::endl;
      exit(1);
    }
#endif
    unsigned outputAgents = parameterMap.at("OUTPUT_AGENTS_DURING_SIM").dbl(0);
    unsigned outputFrequency =
      parameterMap.at("OUTPUT_AGENTS_DURING_SIM").dbl(1);
    bool analyzeAfterStabilization =
      parameterMap.at("ANALYZE_AFTER_STABILIZATION").isSet();
    unsigned analyzeFrequency = parameterMap.at("ANALYZE_DURING_SIM").dbl();
    unsigned outputAgentsAfterStabilization =
      parameterMap.at("OUTPUT_AGENTS_AFTER_STABILIZATION").dbl();
    for (unsigned i = 0; i < numIterations; ++i, currentDate += timeStep) {
      for (auto& e: events) e(this);
      if (timing > 0 && (i + 1) % timing == 0) {
        gettimeofday(&timeEnd, NULL);
        elapsedTime = timeEnd.tv_sec - timeBegin.tv_sec;
        csvout("TIMING","DURING", elapsedTime);
      }
      if (outputAgents > 0 && (i + 1) % outputFrequency == 0) {
        printAgents(agents, simulationNum, currentDate);
      }
      if ( stabilizationSteps > 0 && stabilizationSteps == i) {
        if (analyzeAfterStabilization) analysis(true, true, true, analyzeTruncatedAgeInit);
        if (outputAgentsAfterStabilization) {
          printAgents(agents, simulationNum, currentDate);
        }
      } else if (analyzeFrequency &&
                 (i + 1) % analyzeFrequency == 0) {
        analysis(false, false, false, analyzeTruncatedAgeDuring);
      }
    }

    gettimeofday(&timeEnd, NULL);
    elapsedTime = timeEnd.tv_sec - timeBegin.tv_sec;
    if (parameterMap.at("OUTPUT_TIMING_AFTER").dbl()) {
      csvout("TIMING", "AFTER", elapsedTime);
    }

    /* Wrap up */
    outputAgents = parameterMap.at("OUTPUT_AGENTS_AFTER").dbl();
    if (outputAgents) printAgents(agents, simulationNum, endDate);
    if ( (unsigned) parameterMap.at("ANALYZE_AFTER").isSet()) {
      analysis(true, true, true);
    }
  }

  /**
     Populates a matrix of doubles from a CSV file.

     @param key[in] ParameterMap key that contains name of CSV file
     @param delim[in] Delimiter (usually semi-colon or comma) of the csv file
     @param hasHeader[in] Whether or not the CSV file has a header
   */

  inline DblMatrix matrixFromCSV(const char* key,
                                 const char* delim,
                                 const bool hasHeader)
  {
    std::string filename = parameterMap.at(key).str();
    CSVParser csvParser(filename.c_str(), delim, hasHeader);
    DblMatrix result = csvParser.toDoubles();
    return result;
  }

  /**
      Creates agents using Stefan Scholz's algorithm. Ported from R, hence
      some peculiarities. This function needs to be neatened and optimized.
   */
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

    size_t rows = maxAge - MIN_AGE;
    truncateMatrix(demographics,maxAge, 7);
    truncateMatrix(singles, rows, 5);
    truncateMatrix(partners, rows, 5);

    assert(singles.size());

    // Calculate initial infection rates matrices
    DblMatrix initialRates = matrixFromCSV("INITIAL_INFECTION_RATES_CSV",
                                           ",", true);
    double scaleInitialRates = parameterMap.at("SCALE_INITIAL_RATES").dbl();
    std::vector<double> initialInfectionRatesMSW(MAX_AGE, 0.0);
    std::vector<double> initialInfectionRatesMSM(MAX_AGE, 0.0);
    std::vector<double> initialInfectionRatesWSM(MAX_AGE, 0.0);
    std::vector<double> initialInfectionRatesWSW(MAX_AGE, 0.0);
    {
      unsigned i = 0;
      for (auto& row : initialRates) {
        for (;i <= row[0] && i < MAX_AGE; ++i) {
          initialInfectionRatesMSW[i] = row[1] * scaleInitialRates;
          initialInfectionRatesMSM[i] = row[2] * scaleInitialRates;
          initialInfectionRatesWSM[i] = row[3] * scaleInitialRates;
          initialInfectionRatesWSW[i] = row[4] * scaleInitialRates;
        }
      }
      for (; i < MAX_AGE; ++i) {
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
    // std::vector<double> ageRange;
    // for (unsigned i = 12; i <= 100; ++i) ageRange.push_back(i);
    auto ageShare = getCol(singles, 1);
    auto femRatio = getCol(singles, 2);
    auto msmRate = getCol(singles, 3);
    auto wswRate = getCol(singles, 4);

    // Create the single agents
    createAgents(agents,
                 0, S,
                 ageShare, femRatio,
                 wswRate, msmRate, ww, mw, wm, mm,
                 initialInfectionRatesMSW, initialInfectionRatesMSM,
                 initialInfectionRatesWSM, initialInfectionRatesWSW, false);

    // Get age structure, sex & sexual orientation information for paired agents
    // ageShare = femRatio = msmRate = wswRate = {};
    ageShare = getCol(partners, 1);
    femRatio = getCol(partners, 2);
    msmRate = getCol(partners, 3);
    wswRate = getCol(partners, 4);

    // Create the paired agents and their partners
    createAgents(agents,
                 S, X,
                 ageShare, femRatio,
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


  /**
     Initializes the characteristics of an individual agent.
     Largely ported from Stefan Scholz's algorithm in R, so there
     are pecularities in this code. It needs to be tidied and perhaps
     optimised.
   */
  inline void initAgent(Agent *agent,
                        const unsigned id,
                        Sample<>& sampleAgeshare,
                        const std::vector<double>& femRatio,
                        const std::vector<double>& wswRate,
                        const std::vector<double>& msmRate,
                        const std::vector<double>& initialInfectionRatesMSW,
                        const std::vector<double>& initialInfectionRatesMSM,
                        const std::vector<double>& initialInfectionRatesWSM,
                        const std::vector<double>& initialInfectionRatesWSW,
                        std::vector<Sample<>>& sample_matWW,
                        std::vector<Sample<>>& sample_matMW,
                        std::vector<Sample<>>& sample_matWM,
                        std::vector<Sample<>>& sample_matMM)
  {
    std::uniform_real_distribution<double> uni;

    agent->id = id;

    // Age
    unsigned age = sampleAgeshare() + 12;
    // Add a random part of the year onto the age
    agent->age = age + uni(rng) * YEAR;
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


  /**
     Stochastically infects the partner of an infected agent.
     At the initialization of the simulation, partners of infected agents are
     more likely to be infected than the general population. This method
     implements this idea.

     @param agent[in,out] Partner agent who may be infected.
   */
  inline void conditionalInfectPartner(Agent *agent)
  {
    std::uniform_real_distribution<double> uni;
    if (uni(rng) < probInfectedIfPartnerInfected) agent->infected = true;
  }

  /**
     Creates agents at initialization using code ported from R written
     by Stefan Scholz. Some peculiarities because of the direct "translation"
     remain. This code needs to be cleaned and optimized.
   */

  void createAgents(AgentVector& agents,
                    const unsigned fromAgent,
                    const unsigned toAgent,
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
    Sample<> sample_ageshare(ageShare, &rng);
    vector<Sample<>> sample_matWW(matWW[0].size());
    vector<Sample<>> sample_matMW(matMW[0].size());
    vector<Sample<>> sample_matWM(matWM[0].size());
    vector<Sample<>> sample_matMM(matMM[0].size());
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
        agent->setSinglePeriod(currentDate,
                               weibullSinglePeriodInitial,
                               scaleSinglePeriodInitial,
                               probZeroSinglePeriod,
                               scaleSinglePeriodZeroDaysInitial);
      }
      agents.push_back(agent);
      if (agent->partner) agents.push_back(agent->partner);
    }
  }

  /**
     Calculates population demographics at beginning of simulation, such as
     number of males, females, number infected, number by sexual orientation.
   */
  void calculateDemographics()
  {
    for (unsigned i = 0; i < agents.size(); ++i) {
      if (agents[i]->sex == MALE) {
        ++numMales;
        if (agents[i]->sexual_orientation == HETEROSEXUAL)
          ++numMsw;
        else
          ++numMsm;
        if (agents[i]->infected) {
          ++numInfectedMales;
          if (agents[i]->sexual_orientation == HETEROSEXUAL)
            ++numInfectedMsw;
          else
            ++numInfectedMsm;
        }
      } else {
        ++numFemales;
        if (agents[i]->sexual_orientation == HETEROSEXUAL)
          ++numWsm;
        else
          ++numWsw;
        if (agents[i]->infected) {
          ++numInfectedFemales;
          if (agents[i]->sexual_orientation == HETEROSEXUAL)
            ++numInfectedWsm;
          else
            ++numInfectedWsw;
        }
      }
    }
  }

  /**
     Updates demographics when agent becomes infected.

     @param agent[in] Agent who has become infected.
   */
  void trackRiskFactors(const Agent* agent)
  {
    if (agent->sex == MALE) {
      ++numInfectedMales;
      if (agent->sexual_orientation == HETEROSEXUAL)
        ++numInfectedMsw;
      else
        ++numInfectedMsm;
    } else {
      ++numInfectedFemales;
      if (agent->sexual_orientation == HETEROSEXUAL)
        ++numInfectedWsm;
      else
        ++numInfectedWsw;
    }
  }


  /**
     Creates the mating pool on each iteration of the time step.
   */
  AgentVector getUnmatchedAgents()
  {
    AgentVector matingPool;
    for (auto& agent : agents) {
      if (agent->isMatchable(currentDate))
        matingPool.push_back(agent);
    }
    return matingPool;
  }

  /**
     Shuffles the mating pool.
   */
  AgentVector getShuffledUnmatchedAgents()
  {
    AgentVector matingPool = getUnmatchedAgents();
    shuffle(matingPool.begin(), matingPool.end(), rng);
    // Remove back agent if odd
    if (matingPool.size() % 2 == 1) matingPool.pop_back();
    if (printNumMatings) {
      csvout("MATINGPOOL", "", matingPool.size());
    }
    return matingPool;
  }

  /**
     Calculates the closest partnership in a contiguous list of agents to
     the first one in the list.

     @param from[in] Iterator to agent against which the other agents
     are compared for partnership compatibility.
     @param to[in] Iterator to one past last agent in contiguous list
  */

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

  /**
     Outputs statistics for current date.

     @param orientationStats[in] Whether or not to do sexual orientation stats
     @param ageStats[in] Whether or not to do age structure stats
     @param scoreStats[in] Whether or not to do partnership quality stats
  */

  void analysis(const bool orientationStats = false,
                const bool ageStats = false,
                const bool scoreStats = false,
                const bool truncatedAge = false)
  {
    double prevalence = (double) (numInfectedMales + numInfectedFemales) /
      (numMales + numFemales);
    double malePrevalence = (double) numInfectedMales / numMales;
    double femalePrevalence = (double) numInfectedFemales / numFemales;
    double msmPrevalence = (double) numInfectedMsm / numMsm;
    double wswPrevalence = (double) numInfectedWsw / numWsw;

    // Do whenever we call analysis
    csvout("ANALYSIS", "INFECTED", numInfectedMales + numInfectedFemales);
    csvout("ANALYSIS", "PREVALENCE", prevalence);
    csvout("ANALYSIS", "PARTNERSHIPS", totalPartnerships);

    if (orientationStats) {
      csvout("ANALYSIS", "MALE_PREVALENCE", malePrevalence);
      csvout("ANALYSIS", "FEMALE_PREVALENCE", femalePrevalence);
      csvout("ANALYSIS", "MSM_PREVALENCE", msmPrevalence);
      csvout("ANALYSIS", "WSW_PREVALENCE", wswPrevalence);
      csvout("ANALYSIS", "MSM_PARTNERSHIPS", totalMsmPartnerships);
      csvout("ANALYSIS", "WSW_PARTNERSHIPS", totalWswPartnerships);
    }

    if (ageStats) {
      struct AgeStatistics a;
      calculateAgeStatistics(a);

      for (unsigned i = 0; i < NUM_INTERVALS; ++i) {
        unsigned from = i * AGE_INTERVAL;
        unsigned to = i * AGE_INTERVAL + AGE_INTERVAL - 1;
        std::stringstream stream;
        stream << "_" << from << "-" << to;
        csvout("ANALYSIS", std::string("MALE_AGE") + stream.str(),
               a.malesByAge[i]);
        csvout("ANALYSIS", std::string("FEMALE_AGE") + stream.str(),
               a.femalesByAge[i]);
        ostringstream ssmale, ssfemale;
        if (a.malesByAge[i] == 0)
          ssmale << "NA";
        else
          ssmale << std::setprecision(3)
                 << ( (double) a.infectedMalesByAge[i] / a.malesByAge[i] );
        if (a.femalesByAge[i] == 0)
          ssfemale << "NA";
        else
          ssfemale << std::setprecision(3)
                   << ((double) a.infectedFemalesByAge[i] / a.femalesByAge[i]);
        csvout("ANALYSIS",
               std::string("MALE_PREVALENCE_AGE") + stream.str(), ssmale);
        csvout("ANALYSIS",
               std::string("FEMALE_PREVALENCE_AGE") + stream.str(), ssfemale);
      }
    }

    if (scoreStats) {
      csvout("ANALYSIS", "SCORE", totalPartnershipScore / totalPartnerships);
      csvout("ANALYSIS", "FAILED", failedMatches);
      csvout("ANALYSIS", "POOR", poorMatches);
    }
    if (truncatedAge) {
      unsigned numAgents = 0;
      unsigned numMaleAgents = 0;
      unsigned numFemaleAgents = 0;
      unsigned numMswAgents = 0;
      unsigned numWsmAgents = 0;
      unsigned numMsmAgents = 0;
      unsigned numWswAgents = 0;

      unsigned numInfectedAgents = 0;
      unsigned numInfectedMaleAgents = 0;
      unsigned numInfectedFemaleAgents = 0;
      unsigned numInfectedMswAgents = 0;
      unsigned numInfectedWsmAgents = 0;
      unsigned numInfectedMsmAgents = 0;
      unsigned numInfectedWswAgents = 0;
      double fromAge = minAge;
      double toAge = maxAge;
      if (incTruncatedAge == true) {
        if ((currentDate + EPSILON)  >= startDate) {
          double dateIncrement = currentDate - startDate;
          fromAge += dateIncrement;
          toAge += dateIncrement;
        }
      }
      for (auto& a: agents) {
        if ( (a->age + EPSILON) >= fromAge && a->age < toAge) {
          ++numAgents;
          if (a->sex == MALE) {
            ++numMaleAgents;
            if (a->sexual_orientation == HETEROSEXUAL) {
              ++numMswAgents;
            } else {
              ++numMsmAgents;
            }
          } else {
            ++numFemaleAgents;
            if (a->sexual_orientation == HETEROSEXUAL) {
              ++numWsmAgents;
            } else {
              ++numWswAgents;
            }
          }
          if (a->infected) {
            ++numInfectedAgents;
            if (a->sex == MALE) {
              ++numInfectedMaleAgents;
              if (a->sexual_orientation == HETEROSEXUAL) {
                ++numInfectedMswAgents;
              } else {
                ++numInfectedMsmAgents;
              }
            } else {
              ++numInfectedFemaleAgents;
              if (a->sexual_orientation == HETEROSEXUAL) {
                ++numInfectedWsmAgents;
              } else {
                ++numInfectedWswAgents;
              }
            }
          }
        }
      }
      csvout("ANALYSIS", "AGE_RANGE", numAgents);
      csvout("ANALYSIS", "AGE_RANGE_MALE", numMaleAgents);
      csvout("ANALYSIS", "AGE_RANGE_FEMALE", numFemaleAgents);
      csvout("ANALYSIS", "AGE_RANGE_MSW", numMswAgents);
      csvout("ANALYSIS", "AGE_RANGE_WSM", numWsmAgents);
      csvout("ANALYSIS", "AGE_RANGE_MSM", numMsmAgents);
      csvout("ANALYSIS", "AGE_RANGE_WSW", numWswAgents);
      csvout("ANALYSIS", "AGE_RANGE_PREVALENCE",
             (double) numInfectedAgents / numAgents);
      csvout("ANALYSIS", "AGE_RANGE_PREVALENCE_MALE",
             (double) numInfectedMaleAgents / numMaleAgents);
      csvout("ANALYSIS", "AGE_RANGE_PREVALENCE_FEMALE",
             (double) numInfectedFemaleAgents / numFemaleAgents);
      csvout("ANALYSIS", "AGE_RANGE_PREVALENCE_MSW",
             (double) numInfectedMswAgents / numMswAgents);
      csvout("ANALYSIS", "AGE_RANGE_PREVALENCE_WSM",
             (double) numInfectedWsmAgents / numWsmAgents);
      csvout("ANALYSIS", "AGE_RANGE_PREVALENCE_MSM",
             (double) numInfectedMsmAgents / numMsmAgents);
      csvout("ANALYSIS", "AGE_RANGE_PREVALENCE_WSW",
             (double) numInfectedWswAgents / numWswAgents);
    }
  }

  /**
     Calculates statistics for each age interval.

     @param AgeStatistics[out] Struct to hold statistics
   */

  void calculateAgeStatistics(AgeStatistics& a)
  {
    for (unsigned i = 0; i < agents.size(); ++i) {
      unsigned age = agents[i]->age;
      unsigned interval = age / AGE_INTERVAL;
      if (agents[i]->sex == MALE) {
        ++a.malesByAge[interval];
        if (agents[i]->infected) ++a.infectedMalesByAge[interval];
      } else {
        ++a.femalesByAge[interval];
        if (agents[i]->infected) ++a.infectedFemalesByAge[interval];
      }
    }
  }

  /**
     Writes analytical information in comma separated format to stdout.

     @param desc1[in] Description 1
     @param desc2[in] Description 2
     @param extra[in] Additional info to output
   */
  inline void csvout(const std::string& desc1,
                     const std::string& desc2,
                     const std::ostringstream& extra)
  {
    writeCsvLine(simulationName, simulationNum, currentDate,
                 desc1, desc2, extra);
  }

  /**
     Writes analytical information in comma separated format to stdout.

     @param desc1[in] Description 1
     @param desc2[in] Description 2
     @param value[in] Real number (double) to write out
     @param precision[in] Number of decimal places to write
  */
  inline void csvout(const std::string& desc1,
                     const std::string& desc2,
                     const double value,
                     const unsigned precision = 3)
  {
    std::ostringstream stream;
    stream << std::setprecision(precision) << value;
    csvout(desc1, desc2, stream);
  }

  /**
     Writes analytical information in comma separated format to stdout.

     @param desc1[in] Description 1
     @param desc2[in] Description 2
     @param value[in] Unsigned int to write
  */
  inline void csvout(const std::string& desc1,
                     const std::string& desc2,
                     const unsigned value)
  {
    std::ostringstream stream;
    stream << value;
    csvout(desc1, desc2, stream);
  }

  /**
     Writes analytical information in comma separated format to stdout.

     @param desc1[in] Description 1
     @param desc2[in] Description 2
     @param value[in] C++ size_t to write
  */
  inline void csvout(const std::string& desc1,
                     const std::string& desc2,
                     const size_t value)
  {
    std::ostringstream stream;
    stream << value;
    csvout(desc1, desc2, stream);
  }

  /**
     Makes two agents partners.

     @param a: One of two agents to make partners.
     @param b: Other of two agents to make partners.
     @param score: Compatibility of agents (used to update stats)
   */

  void makePartner(Agent* a, Agent *b,
                   const double score)
  {
    //// assert(a->partner == NULL);
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

  /**
     Measures compatibility of two agents for matching.

     @param a: One of two agents to measure distance.
     @param b: Other of two agents to measure distance.
   */
  double distance(const Agent *a, const Agent *b) const
  {
    return distanceMethod == HEURISTIC_DISTANCE
      ? heuristicDistance(a, b) : tableDistance(a, b);
  }

  /**
     One of the methods for measuring compatibility of two agents for matching.
     This one relies on agent's specific age preference.

     @param a: One of two agents to measure distance.
     @param b: Other of two agents to measure distance.
   */


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

  /**
     One of the methods for measuring compatibility of two agents for matching.
     This one looks up agent's age preference in a table.

     @param a: One of two agents to measure distance.
     @param b: Other of two agents to measure distance.
   */


  double tableDistance(const Agent *a, const Agent *b) const
  {
    double score = 0.0;
    unsigned a_age = std::min( (unsigned) a->age, (unsigned) MAX_AGE - 1) - MIN_AGE;
    unsigned b_age = std::min( (unsigned) b->age,  (unsigned) MAX_AGE - 1) - MIN_AGE;

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

  /**
     Calculates the cluster value for an agent. Useful in some pair-matching
     algorithms, such as CSPM.

     @param a[in] Agent whose cluster value to calculate
   */
  double clusterValue(const Agent *a) const
  {
    return ( (double) a->desired_age / 100.0 + a->age / 100.0) / 2.0
      + a->sexual_orientation;
  }

  virtual void setEvents();

};

#endif
