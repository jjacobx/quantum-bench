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
#define DEPTH_DEFAULT 16

using namespace std;

void random_circuit(Qureg qureg, int depth);

void print_qureg(Qureg qureg);
void set_args(int argc, char *argv[], int &qreg_size, int &depth);
int set_verbose();


int main(int argc, char *argv[]) {
  // Set number of qubits and verbosity
  int qreg_size = QREG_DEFAULT;
  int depth = DEPTH_DEFAULT;

  set_args(argc, argv, qreg_size, depth);
  int verbose = set_verbose();

  // Prepare the hardware-agnostic QuEST environment
  QuESTEnv env = createQuESTEnv();

  if (env.rank == 0 && verbose) {
    cout << "Verbose is ON" << endl;
    cout << "No. processes: " << env.numRanks << endl; 
    cout << "No. qubits: " << qreg_size << endl;
    cout << "Depth: " << depth << endl;
  }

  Qureg qureg = createQureg(qreg_size, env);
  initZeroState(qureg);

  // Sync and start timer
  syncQuESTEnv(env);
  auto tstart = chrono::steady_clock::now();

  // Run QFT
  random_circuit(qureg, depth);
  
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
  
  // Free memory
  destroyQureg(qureg, env);
  destroyQuESTEnv(env);

  return 0;
}


void crot(Qureg qureg, int targetQubit, int controlQubit, int k) {
  unsigned long frac = 1 << k;
  controlledPhaseShift(qureg, targetQubit, controlQubit, 2 * M_PI / frac);
}

void rand(Qureg qureg, int qubit) {
  int r = rand() % 3;
  Vector v;
    if (r == 0) {
      // X axis
      v.x = 1; v.y = 0; v.z = 0;
    } else if (r == 1) {
      // Y axis
      v.x = 0; v.y = 1; v.z = 0;
    } else {
      // X + Y axis
      v.x = 1; v.y = 1; v.z = 0;
    }

    rotateAroundAxis(qureg, qubit, M_PI / 2, v);
}

void random_circuit(Qureg qureg, int depth) {
  for (int i = 0; i < depth; i++) {
    for (int j = 0; j < qureg.numQubitsRepresented; j++)
      rand(qureg, j);
    for (int j = i % 2; j < qureg.numQubitsRepresented - 1; j += 2)
      crot(qureg, j, j + 1, 2);
  }
}

void print_qureg(Qureg qureg) {
  int numStates = 1 << qureg.numQubitsRepresented;
  for (int i = 0; i < numStates; i++) {
    Complex amp = getAmp(qureg, i);
    bitset<QREG_MAX> ibits(i);

    cout << "Amplitude of |" << ibits << ">: ";
    printf("%.8f + %.8fi\n", amp.real, amp.imag);
  }
}

void set_args(int argc, char* argv[], int& qreg_size, int& depth) {
  for (int i = 1; i < argc; i++) {
    string arg = argv[i];
    if (arg == "-q") {
      qreg_size = atoi(argv[++i]);
    } else if (arg == "-d") {
      depth = atoi(argv[++i]);
    } else {
      string message = "Error: Unknown argument '" + arg + 
        "'! Use: ./bin -q $NQUBITS -d $DEPTH";
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
