#include <iostream>
class Heisenberg
{
private:
 double** m_chain;
 double* m_TrialSpin;
 int m_nodes;
 int max_iter;
 double k_temp;
 double m_j;
 double m_smin;
public:
 Heisenberg(int nodes=40, double temp=1, double j=1, double s_min=0.2, int max_iter=4000);
 ~Heisenberg();
 double Transition(double d_energy);
 double CalculateEnergy();
 double max(double x, double y, double z);
 double deltaS();
 void MonteCarloMove(double s_max);
 double random(double m, double n);
 void setInitialState();
 void trialSpin(int element);
 void printLines();

};
#endif
