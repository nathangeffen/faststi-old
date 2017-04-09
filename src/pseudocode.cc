/* Agent data structure */

class Agent {
  partner; /* Pointer to agent */
  past_partners; /* List of part partners */
  infectivity; /* How easily infected this agent is [0..1] */
  infectiousness; /* How easily this agent infects others [0..1] */
  breakupiness; /* Measure of how easily an agent breaks up with partner */
  infected; /* Bool indicating whether infected */
  age; /* Age of agent 12-100 */
  sex; /* MALE or FEMALE */
  preferred_age; /* Ideal age of partner */
  preferred_sex; /* Ideal sex of partner */
  wants_relationship; /* true or false */
};


/******* Parameters *************/
/* These are user specified. Values below are defaults */

NUM_AGENTS = 10000000;
MALE_INFECTIVITY = 0.1;
FEMALE_INFECTIVITY = 0.2;
MALE_INFECTIOUSNESS = 0.2;
FEMALE_INFECTIOUSNESS = 0.1;
BREAKUPINESS = 0.1; /* This is perhaps way too high */
START_DATE = 2017;
END_DATE = 2018;
TIME_STEP = 1.0 / 365.0;
MATCHING_ALGORITHM = one of [BFPM, RPM, CSPM];
PREV_PARTNER_PENALTY = 10.0;


/******/

/* Engine */

void simulate()
{
  Agent agents[] = initPopulation();
  for (i = START_DATE; i < END_DATE; i = i + TIME_STEP) {
    for (each event in [age, breakup, pair, infect]) {
      event(agents);
    }
  }
}

initPopulation()
{
  Agent *agents;
  /* Stefan: Noting your suggestion in the analysis plan -
     But perhaps we could simply initiate the population with
     the C++ version of your algorithm?
  */
  /* .... After creating the agents and their partnerships */
  /* Initialize the remaining parameters */

  for (agent in agents) {
    /* VERSION 1 */
    if (agent.sex == MALE) {
      agent.infectivity = MALE_INFECTIVITY;
      agent.infectiousness = MALE_INFECTIOUSNESS;
    } else {
      agent.infectivity = FEMALE_INFECTIVITY;
      agent.infectiousness = FEMALE_INFECTIOUSNESS;
    }
    agent.breakupiness = BREAKUPINESS;
  }
  return agents;
}

/* Events */

void age(agents)
{
  for (agent in agents) {
    agent.age = agent.age + TIME_STEP;
  }
}

void breakup(agent)
{
  for (agent in agents) {
    if (agent has partner) {
      risk_breakup = (agent.breakupiness + agent.partner.breakupiness) / 2;
      if (uniform random number in [0..1] < risk_breakup) {
        agent.partner.past_partners.append(agent);
        agent.partner.partner = null;
        agent.past_partners.append(agent.partner);
        agent.partner = null;
      }
    }
  }
}

void pair(agent)
{
  MATCHING_ALGORITHM(agent);
}

void infect(agents)
{
  for (agent in agents) {
    if (agent is not infected and agent has partner and partner is infected) {
      risk_infection = (agent.infectivity + agent.partner.infectiousness) / 2;
      if (uniform random number in [0..1] < risk_infection) {
        agent.infected = true;
      }
    }
  }
}

/* Algorithms */

void BFPM(agents);
void CSPM(agents);
void RSPM(agents);

/* Support functions */

void approx_distance(a, b)
{
  score = 0; /* The lower the score the better. */
  if (a.wants_relationship == false) {
    score += 1000;
  }
  if (b.wants_relationship == false) {
    score += 1000;
  }
  if (a.preferred_sex != b.sex) {
    score += 1000;
  }
  if (b.preferred_sex != a.sex) {
    score += 1000;
  }
  if (a in b.partner_list) {
    score += PREV_PARTNER_PENALTY;
  }
  age_diff = fabs(a.preferred_age - b.age);
  /* Stefan, I'm open to ideas on what this looks like. */
  score += some distribution function around age_diff;
  /* Repeat for b, a */
  age_diff = fabs(b.preferred_age - a.age);
  score += some distribution function around age_diff;
  return score;
}

void exact_distance(a, b)
{
  /* Same as it currently is in germany_partners.cc */
}
