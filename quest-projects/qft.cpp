#include <bitset>
#include <chrono>
#include <cmath>
#include <iostream>

#include "QuEST.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define QREG_DEFAULT 24
#define QREG_MAX 32
#define PRECISION 1E-10

using namespace std;


void qft(Qureg qureg);

bool operator==(Complex const& lhs, Complex const& rhs);
void validate_result(QuESTEnv env, Qureg& qureg);
void print_qureg(Qureg qureg);

void set_args(int argc, char* argv[], int& qreg_size);
int set_verbose();


int main(int argc, char *argv[]) {
  // Set number of qubits and verbosity
  int qreg_size = QREG_DEFAULT;

  set_args(argc, argv, qreg_size);
  int verbose = set_verbose();

  // Prepare the hardware-agnostic QuEST environment
  QuESTEnv env = createQuESTEnv();

  if (env.rank == 0 && verbose) {
    cout << "Verbose is ON" << endl;
    cout << "No. processes: " << env.numRanks << endl; 
    cout << "No. qubits: " << qreg_size << endl;
  }

  Qureg qureg = createQureg(qreg_size, env);
  initZeroState(qureg);

  // Sync and start timer
  syncQuESTEnv(env);
  auto tstart = chrono::steady_clock::now();

  // Run QFT
  qft(qureg);
  
  // Sync and stop timer
  syncQuESTEnv(env);
  auto tstop = chrono::steady_clock::now();

  double tdiff = chrono::duration<double, milli>(tstop - tstart).count();

  if (env.rank == 0) {
    if (verbose) 
      cout << "Time taken: " << tdiff << " ms" << endl;
    else 
      cout << tdiff << endl;
  }

  // Validate whether all coefficients are the same
  if (verbose)
    validate_result(env, qureg);
  
  // Free memory
  destroyQureg(qureg, env);
  destroyQuESTEnv(env);

  return 0;
}


void crot(Qureg qureg, int targetQubit, int controlQubit, int k) {
  unsigned long frac = 1 << k;
  controlledPhaseShift(qureg, targetQubit, controlQubit, 2 * M_PI / frac);
}

void multi_crot(Qureg qureg, int targetQubit) {
  for (int i = 2; i <= qureg.numQubitsRepresented - targetQubit; i++)
    crot(qureg, targetQubit, targetQubit + i - 1, i);
}

void swap_qureg(Qureg qureg) {
  int n = qureg.numQubitsRepresented;
  for (int i = 0; i < n / 2; i++)
    swapGate(qureg, i, n - i - 1);
}

void qft(Qureg qureg) {
  for (int i = 0; i < qureg.numQubitsRepresented; i++) {
    hadamard(qureg, i);
    multi_crot(qureg, i);
  }

  swap_qureg(qureg);
}


bool operator==(Complex const& lhs, Complex const& rhs) {
  return (fabs(lhs.real - rhs.real) < PRECISION) &&
    (fabs(lhs.imag - rhs.imag) < PRECISION);
}

void validate_result(QuESTEnv env, Qureg& qureg) {
  unsigned long numStates = 1 << qureg.numQubitsRepresented;
  Complex ampZero = getAmp(qureg, 0);
  bool isValid = true;

  for (int i = 1; i < numStates && isValid; i++)
    isValid = getAmp(qureg, i) == ampZero;

  if (env.rank == 0) {
    if (isValid) 
      cout << "Result valid" << endl;
    else
      cout << "Result invalid" << endl;
  }
}

void print_qureg(Qureg qureg) {
  int numStates = 1 << qureg.numQubitsRepresented;
  for (int i = 0; i < numStates; i++) {
    Complex amp = getAmp(qureg, i);
    bitset<QREG_MAX> ibits(i);

    cout << "Probability of |" << ibits << ">: ";
    printf("%.8f + %.8fi\n", amp.real, amp.imag);
  }
}


void set_args(int argc, char* argv[], int& qreg_size) {
  for (int i = 1; i < argc; i++) {
    string arg = argv[i];
    if (arg == "-q") {
      qreg_size = atoi(argv[++i]);
    } else {
      string message = "Error: Unknown argument '" + arg + 
        "'! Use: ./bin -q $NQUBITS";
      throw invalid_argument(message);
    }
  }
}

int set_verbose() {
  const char* tmp = getenv("VERBOSE");
  string verbose_str(tmp ? tmp : "");
  if (verbose_str == "1")
    return 1;
  else
    return 0;
}
