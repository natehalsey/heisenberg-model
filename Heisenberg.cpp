#include "Heisenberg.h"
#include <cmath>
#include <fstream>

Heisenberg::Heisenberg(int nodes, double temp, double j, double s_min, int max_iter)
{
    m_nodes = nodes;
    k_temp = temp;
    m_j = j;
    m_smin = s_min;
    
    m_TrialSpin = new double[3];
    m_chain = new double*[m_nodes];
    for (int i = 0; i < m_nodes; i ++)
    {
        m_chain[i] = new double[3];
    }

}
Heisenberg::~Heisenberg()
{
    for (int i = 0; i < m_nodes; i ++)
    {
        delete[] m_chain[i];
    }
    delete[] m_chain;
    delete[] m_TrialSpin;
}

double Heisenberg::Transition(double d_energy)
{
    double value;
    double result;
    value = exp(-d_energy/k_temp);
    result = value/(1+value);
    return result;
}
double Heisenberg::CalculateEnergy()
{
    double sum = 0;
    for (int i = 0; i < m_nodes-1; i ++)
    {
        for (int k = 0; k < 3; k ++)
        {
            sum += m_chain[i][k]*m_chain[i+1][k];  
        }
    }    
    for (int j = 0; j < 3; j++)
    {
        sum += m_chain[m_nodes-1][j]*m_chain[0][j];
    }

    return -0.5*m_j*sum;
}
void Heisenberg::setInitialState()
{


    for (int i = 0; i < m_nodes; i ++)
    {
        double a,b,c;
        a = rand();
        b = rand();
        c = rand();
        double magnitude = 1/(sqrt(pow(a,2)+pow(b,2)+pow(c,2)));
        m_chain[i][0] = a*magnitude;
        m_chain[i][1] = b*magnitude;
        m_chain[i][2] = c*magnitude;
    }
}
   
void Heisenberg::MonteCarloMove(double s_max)
{

    double a = random(-1,1);
    double b = random(-1,1);
    double c = random(-1,1);

    a = s_max*a;
    b = s_max*b;
    c = s_max*c;

    if (sqrt(pow(a,2)+pow(b,2)+pow(c,2))>s_max)
    {
        MonteCarloMove(s_max);
    }
    else
    {
        m_TrialSpin[0] = a;
        m_TrialSpin[1] = b;
        m_TrialSpin[2] = c;
    }
    
}

double Heisenberg::random(double m, double n)
{
    return m + (rand() / ( RAND_MAX / (n-m)));
}
void Heisenberg::trialSpin(int element)
{
    double EnergyOld = CalculateEnergy();



    MonteCarloMove(m_smin);
    double a = m_chain[element][0] + m_TrialSpin[0];
    double b = m_chain[element][1] + m_TrialSpin[1];
    double c = m_chain[element][2] + m_TrialSpin[2];

    double magnitude = 1/sqrt(pow(a,2)+pow(b,2)+pow(c,2));

    a = a * magnitude;
    b = b * magnitude;
    c = c * magnitude;

    double holder_a = m_chain[element][0];
    double holder_b = m_chain[element][1];
    double holder_c = m_chain[element][2];

    m_chain[element][0] = a;
    m_chain[element][1] = b;
    m_chain[element][2] = c;
    

    double EnergyNew = CalculateEnergy();

    double deltaEnerg = EnergyOld - EnergyNew;
    double randomNum = random(0,1);

    if (((deltaEnerg) > 0)&&(Transition(deltaEnerg)<randomNum))
    {
        m_chain[element][0] = holder_a;
        m_chain[element][1] = holder_b;
        m_chain[element][2] = holder_c;
    }
    
}
void Heisenberg::printLines()
{

    setInitialState();
    std::ofstream myFile;
    myFile.open("Project.dat");
    for (int i = 0; i < max_iter; i ++)
    {
        for (int j = 0; j < m_nodes; j++)
        {
            trialSpin(j);
            myFile << i << "\t" << CalculateEnergy() << "\n";
        }
    }
    myFile.close();
}
int main(int argc, const char* argv[]) 
{
    int i = 0;
    double j = 0;
    std::string s = "J";
    while (i < 10)
    {
        
        char str[10];
        sprintf(str,"%d.txt",i);
        std::string filename = s;
        filename.append(str);
        std::ofstream myFile;

        myFile.open(filename.c_str());

        Heisenberg myHeis = Heisenberg(40,1,j,0.15,4000);
        myHeis.setInitialState();
        
        
        for (int i = 0; i < 4000; i ++)
        {
            for (int j = 0; j < 40; j++)
            {
                myHeis.trialSpin(j);
            }
            myFile << i << "\t" << myHeis.CalculateEnergy() << "\n";
        }
        myFile.close();
        i += 1;
        j += 0.1;
    }
}
