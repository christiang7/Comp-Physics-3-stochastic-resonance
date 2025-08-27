/*
Computational Physics Project
stochastical resonance

Alexander Putz 763265
Christian Gößl 762627

*/

#include <iostream>
#include <math.h>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <random>
#include <ctime>
#include <cstdlib>
#include <vector>

using namespace std;

//declaration of used variables

default_random_engine gen;					//random number generators (gauss distribution)
normal_distribution<double> gauss(0, 1.0);	// gauss(mean, standard deviation)

ofstream stream;  //file outputstream
ofstream stream_T_K;  // T_K file outputstream
ofstream stream_a_k; //for writing a_k
ofstream stream_b_k; //for writing b_k
ofstream stream_ampl_D;
ofstream stream_A_phi;

double tstep;	//steps for simulation
double M_PII = 3.14159265358979323846;
double A;
double omega;
double D;
double T;
double L;
double N;

int averages = 10; // std=10		number of averagings for Ampl-D relation
int tempdown = 0;
int tempup = 0;
double ampl_init;
double ampl_end;
double D_init;
double D_end;
double delta_a = 0;
double T_K0; // original Kramers time scale
double T_K; // Kramers time scale calculated with the frequency of the induced oscillation
bool mess = false; // indicator of states up and down
bool switchi = true; // switch for output yes or no
vector<long double> x;	//dynamic x Array
vector<long double> T_Kdown; // counts the down states
vector<long double> T_Kup; // count sthe up states
long double number;
// control of output of datapoints
double n_of_datapoints_to_file = 200;	//number of datapoints that shall be written to file
double GaussNoise(double mu, double q)
{
   static const double epsilon = std::numeric_limits<double>::min();
   static const double two_pi = 2.0*3.14159265358979323846;

   thread_local double z1;
   thread_local bool generate;
   generate = !generate;

   if (!generate)
      return z1 * q + mu;

   double u1, u2;
   do
   {
      u1 = rand() * (1.0 / RAND_MAX);
      u2 = rand() * (1.0 / RAND_MAX);
   } while (u1 <= epsilon);

   double z0;
   z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
   z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
   return z0 * q + mu;
}
void close_all_streams()
{
   stream.close();
   stream_T_K.close();
   stream_a_k.close();
   stream_b_k.close();
   stream_ampl_D.close();
   stream_A_phi.close();
}

void init(int tempL, double tempA, double tempomega, double tempD, int tempN, bool tempswitchi) {
   close_all_streams();
   srand(time(NULL));  // seed the random number generator with time
   A = tempA; //
   omega = tempomega; // 
   D = tempD; //
   N = tempN; //
   L = tempL;
   tstep = double(L * 2 * M_PII / (N * omega));
   switchi = tempswitchi;
   //stream.open("data_D" + to_string(D) + ".txt");
   stream.open("data.txt");
   stream_T_K.open("T_K.txt", ios::app);
   stream_a_k.open("a_k_omega" + to_string(omega) + ".txt", ios::app);
   stream_b_k.open("b_k_omega" + to_string(omega) + ".txt", ios::app);
   stream_A_phi.open("A_k-phi_k_omega" + to_string(omega) + ".txt", ios::app);
   stream_ampl_D.open("ampl_D.txt", ios::app);
   if (stream.is_open()) {
      stream << "# t	x	gauss" << endl;	//setting header lines
   }
   else	cout << "ERROR opening file!" << endl; //setting error
   number = GaussNoise(0, 1);//
}

struct FFT_a_b {
   double a_k;		//cos Fourier components
   double b_k;		//sin Fourier components
   double A_k;     // Amplitude of a_k and b_k
   double phi_k;   // Phase of a_k and b_k
};
FFT_a_b FFT(vector<long double> x)
{
   FFT_a_b coefficients;
   long double help = 0.;

   help = 0.;
   for (int i = 0; i < N; i++) //integrate Fourier integral
   {
      help += x[i] * cos(omega*i*tstep);
   }
   help *= 1.0 / double(N);

   coefficients.a_k = help; //add to the list

   help = 0.;	//reset for b_k integration
   for (int i = 0; i < N; i++)
   {
      help += x[i] * sin(omega*i*tstep);
   }
   help *= 1.0 / double(N);

   coefficients.b_k = help; //add to list
   coefficients.A_k = sqrt(pow(coefficients.a_k, 2) + pow(coefficients.b_k, 2));
   coefficients.phi_k = atan(coefficients.a_k / coefficients.b_k);
   return coefficients;
}
void Output(int file, long double data1, long double data2, long double data3 = 0, long double data4 = 0, long double data5 = 0)	//data 3, 4 and 5 are optional
{
   //  t   x gauss  
   if (file == 0) stream << data1 << "	" << data2 << "	" << data3 << endl;
   //  k   a_k
   if (file == 1) stream_a_k << data1 << "	" << data2 << endl;
   //  k   b_k
   if (file == 2) stream_b_k << data1 << "	" << data2 << endl;
   //  D   A_k   phi_k
   if (file == 3) stream_A_phi << data1 << "	" << data2 << "	 " << data3 << endl;
   //  ampl
   if (file == 4) stream_ampl_D << data1 << "	" << data2 << " " << data3 << endl;

   if (file == 5) stream_A_phi << data1 << "	" << data2 << "	" << data3 << "	" << data4 << "	" << data5 << endl;
}
long double NumCalc(long double x, long double number, int i) //stochastical differential equation
{
   return x + (x - pow(x, 3) + A * cos(omega*double(i)*tstep))*tstep + number*sqrt(2 * D * tstep);
}
void kramer(int k, double mean)
{
   int j = k;
   if ((x.at(j) < -0.99) && (mess == false))
   {
      tempdown = j;
      mess = true;
      if ((tempup > 1) || (mean > 0))
      {
         T_Kup.push_back(double(j - tempup)*tstep);
      }
   }
   if ( (x.at(j) > 0.99) && (mess == true) )
   {
      T_Kdown.push_back(double(j - tempdown)*tstep);
      tempup = j;
      mess = false;
   }
   if (j + 1 == N)
   {
      if (mess == true)
         T_Kdown.push_back(double(j - tempdown)*tstep);
      else
         T_Kup.push_back(double(j - tempup)*tstep);
   }
}
void set_Data(double L, double A, double omega, double D, double N)
{
   //L, A, omega, D, N
   init(L, A, omega, D, N, switchi);
   stream << "# parameter(L, A, omega, D, N, tstep): " << L << ", " << A << ", " << omega << ", " << D_init << ", " << D_end << ", " << N <<  ", " << tstep <<endl;
   Output(0, 0., x.at(0), number);
  
   for (int i = 0; i < N; i++) {
      x.push_back(NumCalc(x.at(i), number, i));	//
      if ((i%int(n_of_datapoints_to_file) == 0) and (switchi==true))//write every "reduced" step to data, write only 100 steps to file
        Output(0, double(i + 1)*tstep, x[i + 1], number);	// write calculated data to file
      number = GaussNoise(0,1);// 
   }
}
void docu()
{
   set_Data(30, 0, 0.5, 0.4, 1E8);
   switchi=true;
   ampl_init = A;
   ampl_end = 0.5;
   D_init = D;
   D_end = 0.5;
}
void Kramer_analysis()
{
   remove("T_K.txt");
   docu();
   stream_T_K << "# parameter(L, A, omega, D_init, D_end, N, tstep): " << L << ", " << A << ", " << omega << ", " << D_init << ", " << D_end << ", " << N <<  ", " << tstep <<endl;
   for (double D = D_init; D <= D_end; D = D + 0.02)
   {
      double mean_T_K_size = 0;
      int M = 1;
      for (int m = 1; m <= M; m++) // how many times of repititions
      {
         double mean = 0;
         for (int i = 0; i < N; i++)
         {
            mean = x[i] + mean;
            kramer(i, mean);
         }
         mean_T_K_size = mean_T_K_size + (T_Kup.size() + T_Kdown.size());
         x.clear();
         T_Kup.clear();
         T_Kdown.clear();
         x.push_back(0);
         set_Data(L, A, omega, D, N);
      }
      stream_T_K << (1/D) << " " << log(N*tstep / mean_T_K_size*M) << endl;
   }
}
double return_variance(vector <double> data, double avg) {
   double deviation = 0;

   for (int i = 0; i < data.size(); i++) {	//calculate stochastic variance
      deviation += pow((data[i] - avg), 2);
   }

   return deviation / data.size();

}
void Fourier_task()
{
   double amp_avg = 0;
   double phi_avg = 0;
   vector <double> amp_vec;
   vector <double> phi_vec;


   double amp_dev = 0;
   double phi_dev = 0;

   for (double D = 0; D <= 0.25; D += 0.0025)
   {
      cout << " at D=" << D;
      for (int i = 0; i < averages; i++) {	// averages- Iterations to get statistical data

                                    //L, A, omega, D, N
         set_Data(45, 0.1, 0.001, D, 2E6);

         FFT_a_b coefficients = FFT(x);

         amp_avg += coefficients.A_k;		// sum up all amplitudes and phi to get their average value
         phi_avg += coefficients.phi_k;

         amp_vec.push_back(coefficients.A_k);	// store values in a vector/array to be able to calculate variance
         phi_vec.push_back(coefficients.phi_k);

         //for troubleshooting
         if (i == 1) {
            cout << "D= " << D << endl;
            cout << "a_k= " << coefficients.a_k << endl;
            cout << "b_k= " << coefficients.b_k << endl;

         }
         x.clear();
         x.push_back(0);
      }
      amp_avg /= averages;
      cout << "  amp_avg = " << amp_avg<<endl;
      phi_avg /= averages;

      //get deviation out of amp_vec
      amp_dev = return_variance(amp_vec, amp_avg);
      phi_dev = return_variance(phi_vec, phi_avg);

      Output(5, D, amp_avg, phi_avg, amp_dev, phi_dev);
      amp_avg = 0;
      phi_avg = 0;
      amp_dev = 0;
      phi_dev = 0;

   }
}
int main(int argc, char* argv[]) { 
    
   x.push_back(0);

   bool only_FF = false;
   if (only_FF == true)

      Fourier_task();

   else Kramer_analysis();

   close_all_streams();

   return 0;
}
