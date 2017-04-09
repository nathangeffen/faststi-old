/*
  Notes for Stefan

  1. Thanks. This code is very easy to integrate into the simulation. So once
  the bugs are sorted out, this will be very useful.

  2. Check my main() function below to see how I've tested your code.

  3. I removed the id parameter which you use just to seed the random number
  generator. You only need to seed a random number generator once for the
  duration of its life (or not at all if you want to generate the same results
  every time you run the code). In the final C++ code, I'll be converting rand()
  to the C++ Mersenne Twister anyway.

  4. Check the output of rel_part_ecr. The results aren't quite correct
  yet. Definitely not very Weibullish at the moment. :-)

  Questions for Stefan

  # rel_part_ecr

  1. Is 0 FEMALE and 1 MALE?
  0 is FEMALE

  2. Is 0 HOMOSEXUAL and 1 HETEROSEXUAL?
  0 is HETEROSEXUAL

  3. Is rel_part_ecr returning the number of days in the relationship?
  Gives the day in the simulation at which the relationship will end.

  4. For a new relationship, both agents have to end the relationship on
  the same day. As the code stands at the moment, only one agent is considered.
  Do I call the function twice, once with the characteristics of the one agent
  and once with the other, then take the mean of the two?

  5. There is no accounting here for heterogeneity between agents with respect
  to duration length. Should I multiply the value returned by rel_part_ecr by an
  agent's individual propensity to remain in relationships? If you use the
  agent's id for variation, it can only work if this is a calculation of the
  average duration of that agent's relationship at the beginning of the
  simulation. Merely changing the seed won't work then I'm afraid (though I
  think I can come up with a modification to make it work). Also, it opens the
  question, how do we stochastically change relationship lengths during the
  simulation? Or am I misunderstanding?

  rel_single_ecr

  1. Some of my questions above apply here too.

  2. Is bcr duration of the previous relationship? Is the duration of singlehood
  really dependent on that? What if there has been no previous relationship?

  3. Is rnormal just a standard normal function? Is it the PDF or CDF? I've cut
  and pasted a quick and dirty one I found on Stack Overflow below. When I take
  this over to simulate.cc I'll simply use C++'s standard one. I don't
  understand why you call rnormal(ID). I've modified the code to replace ID with
  an agent id, but that's to get the code to compile - I don't understand what
  I've done.

  4. Not sure what ITERS is. I set it to 100.
  Day in simulation.

 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ITERS 100

// WEIBULL-distribution function according to http://www.taygeta.com/random/weibull.html
// WeibullRandomNumbers = scale.*( -log(1-rand(noOfRandomNumbers,1))).^(1/shape);

double rweibull(double wscale, double wshape){
  double zufall;
  zufall = ((double) rand()/ (double) (RAND_MAX) );
  return(wscale*pow((-log(zufall)),(1/wshape)));
}

double rnormal(double value)
{
   return 0.5 * erfc(-value * M_SQRT1_2);
}

// Function returning the duration of a newly formed relationship in dependence of
//  AGE, SEX, SEXOR

int rel_part_ecr(int arg_sex, int arg_age, int arg_sexor){
  // Initialisation of weibull-distribution (see end of script)
  double mu; double sigma;
  // Calculate parameters of weibull
  mu = 12.66171*arg_age - 3.451597*arg_age*arg_age + 0.3707485*arg_age*arg_age*arg_age - 0.02023791*arg_age*arg_age*arg_age*arg_age + 0.0005982764*arg_age*arg_age*arg_age*arg_age*arg_age - 0.00000913783*arg_age*arg_age*arg_age*arg_age*arg_age*arg_age + 0.00000005655356*arg_age*arg_age*arg_age*arg_age*arg_age*arg_age*arg_age - 1.843016*arg_sexor - 0.2841054 - 1.023666*arg_sex + 1.370621*(arg_sexor*arg_sex);
  sigma = -1.768286*arg_age + 0.4386602*arg_age*arg_age - 0.04593379*arg_age*arg_age*arg_age + 0.002485976*arg_age*arg_age*arg_age*arg_age - 0.0000732176*arg_age*arg_age*arg_age*arg_age*arg_age + 0.000001115761*arg_age*arg_age*arg_age*arg_age*arg_age*arg_age - 0.000000006893729*arg_age*arg_age*arg_age*arg_age*arg_age*arg_age*arg_age + 0.2973045*arg_sexor + 0.1596272 + 0.1161272*arg_sex - 0.1804314*(arg_sexor*arg_sex);
  mu = exp(mu);
  sigma = exp(sigma);
  double dummy1 = rweibull(mu, sigma)*360;
  int ecr = (int) dummy1;
  //printf("dummy is %d, My new ECR is: %d \n", dummy2, bcr);
  return(ecr);
}


// Function returning the duration of the time being single after a break-up in dependence of
//  AGE, SEX, SEXOR, DURATION OF RELATIONSHIP, PERSONAL VARIATION (random effect)


int rel_single_ecr(int id, int arg_bcr, int arg_sex, int arg_age, int arg_sexor){

  // Calculate duration of single time
  double befdur1; double befdur2; double befdur3; double befdur4; double befdur5;
  double befdur6;

  befdur1 = (double)(ITERS - arg_bcr)/365;

  if(befdur1 >= 30) {
    befdur1 = 30;
  }

  befdur2 = befdur1*befdur1; befdur3 = befdur2*befdur1; befdur4 = befdur3*befdur1; befdur5 = befdur4*befdur1; befdur6 = befdur5*befdur1;

  double mu; double sigma;

  // Calculate parameters of weibull
  mu = 31.7957032104883 + 0.243575821828405*arg_sex - 0.16357453436092*arg_sexor + sqrt(0.510540557400565)*rnormal(id) - 6.68293057847908*arg_age + 0.54909519430458*(arg_age*arg_age) - 0.021909227574598*(arg_age*arg_age*arg_age) + 0.000427451560082975*(arg_age*arg_age*arg_age*arg_age) - 0.00000326795090171358*(arg_age*arg_age*arg_age*arg_age*arg_age) + 0.0500964076675074*befdur1 - 0.0281480536301023*befdur2 + 0.00372041809520928*befdur3 - 0.000223759143615591*befdur4 + 0.00000674808532236945*befdur5 - 0.0000000835218749639112*befdur6;

  sigma = -6.20642059757102 - 0.028527416744719*arg_sex + 0.0357281740921485*arg_sexor + 1.11378838652714*arg_age - 0.067938761209021*(arg_age*arg_age) + 0.00182536653429785*(arg_age*arg_age*arg_age) - 0.0000224878375337562*(arg_age*arg_age*arg_age*arg_age) + 0.000000107834215112116*(arg_age*arg_age*arg_age*arg_age*arg_age) + 0.223925971786167*befdur1 - 0.0616047635527968*befdur2 + 0.00984539583466668*befdur3 - 0.000800891252438947*befdur4 + 0.000029931584033257*befdur5 - 0.000000403677855052203*befdur6;
  mu = exp(mu);
  sigma = exp(sigma);
  int dummy1 = rweibull(mu, sigma)*360;
  int ecr = (int) dummy1;
  ////printf("dummy is %d, My new ECR is: %d \n", dummy2, bcr);
  return(ecr);
}



int main(int argc, char *argv)
{

  printf("sex\tage\torient\tecr\n");
  for (int s = 0; s <= 1; ++s) {
    for (int a = 12; a <= 100; ++a) {
      for (int o = 0; o <= 1; ++o) {
        int ecr = rel_part_ecr(s, a, o);
        printf("%d\t%d\t%d\t%d\n", s, a, o, ecr);
      }
    }
  }

  printf("id\tbcr\tsex\tage\torient\tecr\n");
  int id = 50;

  for (int b = 100; b < 200; ++b) {
    for (int s = 0; s <= 1; ++s) {
      for (int a = 12; a <= 100; ++a) {
        for (int o = 0; o <= 1; ++o) {
          int ecr = rel_single_ecr(id, b, s, a, o);
          printf("%d\t%d\t%d\t%d\t%d\t%d\n", id, b, s, a, o, ecr);
        }
      }
    }
  }

  return 0;
}
