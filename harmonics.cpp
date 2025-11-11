#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <algorithm>

//=============================================================================
// Complex Number Class
//=============================================================================

class Complex {
public:
  double real;
  double imag;
  
  Complex(double r = 0.0, double i = 0.0) : real(r), imag(i) {}
  
  Complex operator+(const Complex& other) const {
    return Complex(real + other.real, imag + other.imag);
  }
  
  Complex operator-(const Complex& other) const {
    return Complex(real - other.real, imag - other.imag);
  }
  
  Complex operator*(const Complex& other) const {
    return Complex(
      real * other.real - imag * other.imag,
      real * other.imag + imag * other.real
    );
  }
  
  Complex operator*(double scalar) const {
    return Complex(real * scalar, imag * scalar);
  }
  
  Complex operator/(const Complex& other) const {
    double denom = other.real * other.real + other.imag * other.imag;
    if (denom < 1e-30) {
      return Complex(1e30, 0);  // Avoid division by zero
    }
    return Complex(
      (real * other.real + imag * other.imag) / denom,
      (imag * other.real - real * other.imag) / denom
    );
  }
  
  Complex operator/(double scalar) const {
    return Complex(real / scalar, imag / scalar);
  }
  
  double magnitude() const {
    return std::sqrt(real * real + imag * imag);
  }
  
  double phase() const {
    return std::atan2(imag, real);
  }

  // Complex conjugate
  Complex conjugate() const {
    return Complex(real, -imag);
  }

  // Real power: P = Re(V × I*) for phasors (peak values)
  // Average power = P/2
  static double realPower(const Complex& V, const Complex& I) {
    Complex Iconj = I.conjugate();
    Complex S = V * Iconj;  // Complex power
    return S.real / 2.0;     // Average real power
  }

  // Power delivered to resistive load R (assumes RMS voltage)
  double power(double R) const {
    double Vrms = magnitude();
    return Vrms * Vrms / R;
  }
  
  // Parallel combination of impedances
  static Complex parallel(const Complex& z1, const Complex& z2) {
    return (z1 * z2) / (z1 + z2);
  }
};

//=============================================================================
// Component Base Class
//=============================================================================

class Component {
public:
  virtual ~Component() = default;
  
  // Return impedance at angular frequency omega
  virtual Complex impedance(double omega) const = 0;
  
  // Return admittance
  Complex admittance(double omega) const {
    return Complex(1.0, 0.0) / impedance(omega);
  }
  
  // Power dissipation given RMS current and frequency
  virtual double powerDissipation(double irms) const = 0;
};

//=============================================================================
// Capacitor Component
//=============================================================================

class Capacitor : public Component {
private:
  double c;       // Farads
  double esr;     // Ohms
  double rating; // Volts
  
public:
  Capacitor(double c, double esr, double rating = 1000.0)
    : c(c), esr(esr), rating(rating) {}
  
  Complex impedance(double omega) const override {
    // Z = ESR - j/(ωC)
    double Xc = -1.0 / (omega * c);
    return Complex(esr, Xc);
  }
  
  double powerDissipation(double i) const override {
    // P = I²·ESR
    return i * i * esr;
  }
  
  double voltageStress(double i, double omega) const {
    Complex z = impedance(omega);
    return i * z.magnitude();
  }
  
  bool isOverstressed(double i, double omega) const {
    return voltageStress(i, omega) > rating;
  }
  
  double getCapacitance() const { return c; }
  double getESR() const { return esr; }
};

//=============================================================================
// Inductor Component
//=============================================================================

class Inductor : public Component {
private:
  double L;       // Henries
  double dcr;     // Ohms
  double rating;  // Amps
  
public:
  Inductor(double L, double dcr, double rating = 2.0)
    : L(L), dcr(dcr), rating(rating) {}
  
  Complex impedance(double omega) const override {
    // Z = DCR + jωL
    double Xl = omega * L;
    return Complex(dcr, Xl);
  }
  
  double powerDissipation(double i) const override {
    // P = I²·DCR
    return i * i * dcr;
  }
  
  bool isOverstressed(double i) const {
    return i > rating;
  }
  
  double getInductance() const { return L; }
  double getDCR() const { return dcr; }
};

//=============================================================================
// Linear System Solver (Gaussian Elimination)
//=============================================================================

void solveLinearSystem(Complex A[4][4], Complex b[4], Complex x[4], int n) {
  // Create augmented matrix
  Complex aug[4][5];
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      aug[i][j] = A[i][j];
    }
    aug[i][n] = b[i];
  }
  
  // Forward elimination
  for (int i = 0; i < n; ++i) {
    // Find pivot
    int maxRow = i;
    double maxVal = aug[i][i].magnitude();
    for (int k = i + 1; k < n; ++k) {
      double val = aug[k][i].magnitude();
      if (val > maxVal) {
        maxVal = val;
        maxRow = k;
      }
    }
    
    // Swap rows
    if (maxRow != i) {
      for (int k = 0; k <= n; ++k) {
        Complex temp = aug[i][k];
        aug[i][k] = aug[maxRow][k];
        aug[maxRow][k] = temp;
      }
    }
    
    // Make all rows below this one 0 in current column
    for (int k = i + 1; k < n; ++k) {
      Complex factor = aug[k][i] / aug[i][i];
      for (int j = i; j <= n; ++j) {
        aug[k][j] = aug[k][j] - aug[i][j] * factor;
      }
    }
  }
  
  // Back substitution
  for (int i = n - 1; i >= 0; --i) {
    x[i] = aug[i][n];
    for (int j = i + 1; j < n; ++j) {
      x[i] = x[i] - aug[i][j] * x[j];
    }
    x[i] = x[i] / aug[i][i];
  }
}

//=============================================================================
// 5th-Order Chebyshev LPF
//=============================================================================

class ChebyshevLPF {
private:
  std::unique_ptr<Capacitor> c1, c2, c3, c4;
  std::unique_ptr<Inductor> L1, L2, L3;
  double z0 = 200.0; // System impedance

public:
  struct SolutionPoint {
    Complex Vinput;  // Source voltage (before source impedance)
    Complex Vnode1;  // Input node (after source impedance, at C1/L1 junction)
    Complex Vnode2;  // Between L1 and C2 (at C2 junction)
    Complex Vnode3;  // Between L2 and C3 (at C3 junction)
    Complex Voutput; // Output node (after L3, at C4/load junction)

    Complex Ic1, Ic2, Ic3, Ic4;
    Complex Il1, Il2, Il3;

    double psource;      // DC power from source (includes source impedance loss)
    double pinput;       // Real power into filter (after source impedance)
    double pfundamental; // Real power delivered to load
    double pC1, pC2, pC3, pC4;
    double pL1, pL2, pL3;
  };

  ChebyshevLPF(double C1, double C2, double C3, double C4,
               double L1, double L2, double L3,
               double ESRcap, double DCRind)
    : c1(std::make_unique<Capacitor>(C1, ESRcap, 1000.0)),
      c2(std::make_unique<Capacitor>(C2, ESRcap, 1000.0)),
      c3(std::make_unique<Capacitor>(C3, ESRcap, 1000.0)),
      c4(std::make_unique<Capacitor>(C4, ESRcap, 1000.0)),
      L1(std::make_unique<Inductor>(L1, DCRind, 2.0)),
      L2(std::make_unique<Inductor>(L2, DCRind, 2.0)),
      L3(std::make_unique<Inductor>(L3, DCRind, 2.0))
  { }
  
  SolutionPoint solve(double omega, Complex Vin) {
    // Get admittances at this frequency
    Complex yC1 = c1->admittance(omega);
    Complex yC2 = c2->admittance(omega);
    Complex yC3 = c3->admittance(omega);
    Complex yC4 = c4->admittance(omega);
    Complex yL1 = L1->admittance(omega);
    Complex yL2 = L2->admittance(omega);
    Complex yL3 = L3->admittance(omega);
    Complex yload(1.0 / z0, 0.0);
    // H-bridge + transformer acts as low-impedance voltage source
    // Source impedance ≈ transformer DCR + FET Rds_on ≈ 0.5Ω
    Complex ysource(1.0 / 0.5, 0.0);

    // Build 4×4 matrix for 7th-order shunt-first topology:
    // Vsource --[Zsource]-- V0 --[L1]-- V1 --[L2]-- V2 --[L3]-- V3
    //                       |           |           |           |
    //                      C1          C2          C3          C4, Zload
    //                       |           |           |           |
    //                      GND         GND         GND         GND
    //
    // V0 = input node voltage (after source impedance)
    // V1 = between L1 and L2
    // V2 = between L2 and L3
    // V3 = output node voltage (at load) = Vout

    Complex A[4][4] = {Complex(0,0)};
    Complex b[4] = {Complex(0,0)};
    Complex V[4] = {Complex(0,0)};

    // Node 0 equation (input node after source impedance):
    // (ysource + yC1 + yL1)·V0 - yL1·V1 = ysource·Vin
    A[0][0] = ysource + yC1 + yL1;
    A[0][1] = Complex(0,0) - yL1;
    A[0][2] = Complex(0, 0);
    A[0][3] = Complex(0, 0);
    b[0] = ysource * Vin;

    // Node 1 equation (between L1 and L2):
    // -yL1·V0 + (yL1 + yC2 + yL2)·V1 - yL2·V2 = 0
    A[1][0] = Complex(0,0) - yL1;
    A[1][1] = yL1 + yC2 + yL2;
    A[1][2] = Complex(0,0) - yL2;
    A[1][3] = Complex(0, 0);
    b[1] = Complex(0, 0);

    // Node 2 equation (between L2 and L3):
    // -yL2·V1 + (yL2 + yC3 + yL3)·V2 - yL3·V3 = 0
    A[2][0] = Complex(0, 0);
    A[2][1] = Complex(0,0) - yL2;
    A[2][2] = yL2 + yC3 + yL3;
    A[2][3] = Complex(0,0) - yL3;
    b[2] = Complex(0, 0);

    // Node 3 equation (output node):
    // -yL3·V2 + (yL3 + yC4 + yload)·V3 = 0
    A[3][0] = Complex(0, 0);
    A[3][1] = Complex(0, 0);
    A[3][2] = Complex(0,0) - yL3;
    A[3][3] = yL3 + yC4 + yload;
    b[3] = Complex(0, 0);

    // Solve
    solveLinearSystem(A, b, V, 4);

    // Calculate branch currents
    SolutionPoint sol;
    sol.Vinput = Vin;
    sol.Vnode1 = V[0];  // Input node voltage (not same as Vin!)
    sol.Vnode2 = V[1];  // Between L1 and L2
    sol.Vnode3 = V[2];  // Between L2 and L3
    sol.Voutput = V[3]; // Output node voltage

    // Currents through components (for shunt-first topology):
    sol.Ic1 = yC1 * V[0];  // C1 shunt current
    sol.Ic2 = yC2 * V[1];  // C2 shunt current
    sol.Ic3 = yC3 * V[2];  // C3 shunt current
    sol.Ic4 = yC4 * V[3];  // C4 shunt current
    sol.Il1 = yL1 * (V[0] - V[1]);  // L1 series current
    sol.Il2 = yL2 * (V[1] - V[2]);  // L2 series current
    sol.Il3 = yL3 * (V[2] - V[3]);  // L3 series current

    // Convert phasor magnitudes to RMS
    double ic1rms = sol.Ic1.magnitude() / std::sqrt(2.0);
    double ic2rms = sol.Ic2.magnitude() / std::sqrt(2.0);
    double ic3rms = sol.Ic3.magnitude() / std::sqrt(2.0);
    double ic4rms = sol.Ic4.magnitude() / std::sqrt(2.0);
    double il1rms = sol.Il1.magnitude() / std::sqrt(2.0);
    double il2rms = sol.Il2.magnitude() / std::sqrt(2.0);
    double il3rms = sol.Il3.magnitude() / std::sqrt(2.0);

    sol.pC1 = c1->powerDissipation(ic1rms);
    sol.pC2 = c2->powerDissipation(ic2rms);
    sol.pC3 = c3->powerDissipation(ic3rms);
    sol.pC4 = c4->powerDissipation(ic4rms);
    sol.pL1 = L1->powerDissipation(il1rms);
    sol.pL2 = L2->powerDissipation(il2rms);
    sol.pL3 = L3->powerDissipation(il3rms);

    // Calculate power from DC source (ham radio convention: DC input power)
    // This includes power to filter + power dissipated in source impedance
    Complex Iin = ysource * (Vin - V[0]);  // Input current from source
    sol.psource = Complex::realPower(Vin, Iin);  // Total power from source

    // Calculate real power delivered to the filter input (at node V[0])
    // This is the power that actually enters the filter, after source impedance loss
    sol.pinput = Complex::realPower(V[0], Iin);  // Power at filter input node

    // Calculate real power delivered to load
    Complex Iload = yload * V[3];  // Load current (now at V3)
    sol.pfundamental = Complex::realPower(V[3], Iload);

    return sol;
  }
  
  const Capacitor* getC1() const { return c1.get(); }
  const Capacitor* getC2() const { return c2.get(); }
  const Capacitor* getC3() const { return c3.get(); }
  const Capacitor* getC4() const { return c4.get(); }
  const Inductor* getL1() const { return L1.get(); }
  const Inductor* getL2() const { return L2.get(); }
  const Inductor* getL3() const { return L3.get(); }
};

//=============================================================================
// Square Wave Source with Rise/Fall Time
//=============================================================================

class SquareWaveSource {
private:
  double supply = 60.0;   // Volts (before transformer)
  double tRise = 3e-9;      // 3ns rise time
  double tFall = 3e-9;      // 3ns fall time
  double tRatio = 7.0/3.0;  // Transformer voltage ratio (3:7 turns = 36:196Ω impedance)
  
public:
  SquareWaveSource(double supply = 60.0, double tRise = 3e-9, double tFall = 3e-9)
    : supply(supply), tRise(tRise), tFall(tFall) {}
  
  // Fourier coefficient for nth harmonic of square wave
  // Including finite rise/fall time attenuation
  Complex fourierCoefficient(int n, double f0) {
    if (n == 0) return Complex(0, 0);
    
    // H-bridge differential output cancels even harmonics
    if (n % 2 == 0) return Complex(0, 0);
    
    // Ideal square wave: an = (4V/πn) for odd n
    double amplitude = (4.0 * supply) / (M_PI * n);
    
    // Rise/fall time reduces high frequency content
    // Transfer function: H(ω) = sinc(ω·tRise/2)
    double omegaN = 2.0 * M_PI * n * f0;
    double sincArg = omegaN * tRise / 2.0;
    double reduction = (std::abs(sincArg) > 1e-9) 
      ? std::sin(sincArg) / sincArg 
      : 1.0;
    
    amplitude *= reduction;
    
    // After step-up transformer
    amplitude *= tRatio;
    
    // Return as phasor (peak amplitude, zero phase)
    return Complex(amplitude, 0);
  }
};

//=============================================================================
// Harmonic Analyzer
//=============================================================================

class HarmonicAnalyzer {
public:
  struct ComponentStress {
    double Vpeak;
    double Irms;
    bool overstressed;
    
    ComponentStress() : Vpeak(0), Irms(0), overstressed(false) {}
  };
  
  struct BandResult {
    std::string bandName;
    std::string cornerName;
    double frequency;

    std::vector<double> harmW;	 // Watts at each harmonic
    std::vector<double> harmdBc; // dBc relative to fundamental

    double fundW;
    double totalHarmW;
    double insLossdB;
    double totalFiltDissW;

    // Fundamental frequency power metrics
    double vsourcePeak;      // Source voltage (peak)
    double dcInputPower;     // DC power from Veer
    double rfToFilter;       // RF power into filter
    double fundDissipation;  // Filter dissipation at fundamental

    ComponentStress C1, C2, C3, C4, L1, L2, L3;
    
    // Constructor to initialize vectors
    BandResult(int maxHarm = 10) {
      harmW.resize(maxHarm + 1, 0.0);
      harmdBc.resize(maxHarm + 1, -999.0);
    }
  };
  
  static BandResult analyze(
    ChebyshevLPF& filter,
    SquareWaveSource& source,
    const std::string& bandName,
    const std::string& cornerName,
    double f0,
    int maxHarm = 10)
  {
    BandResult result(maxHarm);
    result.bandName = bandName;
    result.cornerName = cornerName;
    result.frequency = f0;

    double totalDiss = 0.0;
    double totalSourcePower = 0.0;  // DC input power (ham radio convention)
    double totalInPower = 0.0;      // RF power to filter

    // Accumulate stress across all harmonics
    std::vector<double> ic1comps, ic2comps, ic3comps, ic4comps;
    std::vector<double> il1comps, il2comps, il3comps;
    std::vector<double> vc1comps, vc2comps, vc3comps, vc4comps;
    
    // Analyze each harmonic
    for (int n = 1; n <= maxHarm; ++n) {
      double fn = n * f0;
      double omegaN = 2.0 * M_PI * fn;
      
      // Get source voltage at this harmonic
      Complex Vsource = source.fourierCoefficient(n, f0);

      if (Vsource.magnitude() < 1e-6) {
        continue; // Negligible
      }

      // Solve filter at this frequency
      auto sol = filter.solve(omegaN, Vsource);

      // Store fundamental frequency metrics for output
      if (n == 1) {
        result.vsourcePeak = Vsource.magnitude();
        result.dcInputPower = sol.psource;
        result.rfToFilter = sol.pinput;
        result.fundDissipation = sol.pC1 + sol.pC2 + sol.pC3 + sol.pC4 +
                                  sol.pL1 + sol.pL2 + sol.pL3;
      }

      // Accumulate DC input power (ham radio convention)
      totalSourcePower += sol.psource;

      // Accumulate RF power delivered to filter
      totalInPower += sol.pinput;

      // Store output power
      result.harmW[n] = sol.pfundamental;

      // Accumulate component dissipation
      totalDiss += sol.pC1 + sol.pC2 + sol.pC3 + sol.pC4 +
                   sol.pL1 + sol.pL2 + sol.pL3;

      // Accumulate stress (RMS sum for currents and voltages)
      double ic1 = sol.Ic1.magnitude() / std::sqrt(2.0);
      double ic2 = sol.Ic2.magnitude() / std::sqrt(2.0);
      double ic3 = sol.Ic3.magnitude() / std::sqrt(2.0);
      double ic4 = sol.Ic4.magnitude() / std::sqrt(2.0);
      double il1 = sol.Il1.magnitude() / std::sqrt(2.0);
      double il2 = sol.Il2.magnitude() / std::sqrt(2.0);
      double il3 = sol.Il3.magnitude() / std::sqrt(2.0);

      ic1comps.push_back(ic1 * ic1);
      ic2comps.push_back(ic2 * ic2);
      ic3comps.push_back(ic3 * ic3);
      ic4comps.push_back(ic4 * ic4);
      il1comps.push_back(il1 * il1);
      il2comps.push_back(il2 * il2);
      il3comps.push_back(il3 * il3);

      // Voltage stress (peak values)
      vc1comps.push_back(sol.Vnode1.magnitude());
      vc2comps.push_back(sol.Vnode2.magnitude());
      vc3comps.push_back(sol.Vnode3.magnitude());
      vc4comps.push_back(sol.Voutput.magnitude());
    }
    
    // Calculate total RMS currents (RSS sum)
    auto sumRMS = [](const std::vector<double>& components) {
      double sum = 0.0;
      for (double c : components) sum += c;
      return std::sqrt(sum);
    };

    result.C1.Irms = sumRMS(ic1comps);
    result.C2.Irms = sumRMS(ic2comps);
    result.C3.Irms = sumRMS(ic3comps);
    result.C4.Irms = sumRMS(ic4comps);
    result.L1.Irms = sumRMS(il1comps);
    result.L2.Irms = sumRMS(il2comps);
    result.L3.Irms = sumRMS(il3comps);

    // Peak voltages (max across all harmonics)
    result.C1.Vpeak = *std::max_element(vc1comps.begin(), vc1comps.end());
    result.C2.Vpeak = *std::max_element(vc2comps.begin(), vc2comps.end());
    result.C3.Vpeak = *std::max_element(vc3comps.begin(), vc3comps.end());
    result.C4.Vpeak = *std::max_element(vc4comps.begin(), vc4comps.end());

    // Check overstress
    result.C1.overstressed = (result.C1.Vpeak > 1000.0);
    result.C2.overstressed = (result.C2.Vpeak > 1000.0);
    result.C3.overstressed = (result.C3.Vpeak > 1000.0);
    result.C4.overstressed = (result.C4.Vpeak > 1000.0);
    result.L1.overstressed = (result.L1.Irms > 2.0);
    result.L2.overstressed = (result.L2.Irms > 2.0);
    result.L3.overstressed = (result.L3.Irms > 2.0);
    
    // Calculate dBc values
    result.fundW = result.harmW[1];
    
    for (int n = 2; n <= maxHarm; ++n) {
      //      if (result.harmW[n] > 1e-12) {
        result.harmdBc[n] = 10.0 * std::log10(
          result.harmW[n] / result.fundW
        );
	//      }
    }
    
    // Total harmonic power (exclude fundamental)
    result.totalHarmW = 0.0;
    for (int n = 2; n <= maxHarm; ++n) {
      result.totalHarmW += result.harmW[n];
    }
    
    // Insertion loss (using DC input power - ham radio convention)
    if (totalSourcePower > 0) {
      double totalOutPower = 0.0;
      for (int n = 1; n <= maxHarm; ++n) {
        totalOutPower += result.harmW[n];
      }
      result.insLossdB = 10.0 * std::log10(
        totalOutPower / totalSourcePower
      );
    } else {
      result.insLossdB = 0.0;
    }
    
    result.totalFiltDissW = totalDiss;
    
    return result;
  }
};

//=============================================================================
// Band Definition Structure
//=============================================================================

struct BandFilter {
  std::string name;
  double fLow;
  double fHigh;
  double c1;  // Farads
  double l1;  // Henries
  double c2;  // Farads
  double l2;  // Henries
  double c3;  // Farads
  double l3;  // Henries
  double c4;  // Farads

  BandFilter(std::string name,
	     double fLMHz,
	     double fHMHz,
	     double c1pF,
	     double c2pF,
	     double c3pF,
	     double c4pF,
	     double l1nH,
	     double l2nH,
	     double l3nH)
    : name(name),
      fLow(fLMHz*1.0e6),
      fHigh(fHMHz*1.0e6),
      c1(c1pF*1.0e-12),
      l1(l1nH*1.0e-9),
      c2(c2pF*1.0e-12),
      l2(l2nH*1.0e-9),
      c3(c3pF*1.0e-12),
      l3(l3nH*1.0e-9),
      c4(c4pF*1.0e-12)
  { }
};

//=============================================================================
// Tolerance Corner Analysis
//=============================================================================

struct ToleranceCorner {
  std::string name;
  double Cfactor;    // Capacitance multiplier
  double Lfactor;    // Inductance multiplier
  double ESRfactor;  // ESR multiplier
  double DCRfactor;  // DCR multiplier
};

//=============================================================================
// CSV Export Functions
//=============================================================================

void exportResultsToCSV(
  const std::vector<HarmonicAnalyzer::BandResult>& results,
  const std::string& filename)
{
  std::ofstream csv(filename);
  
  if (!csv.is_open()) {
    std::cerr << "Error: Could not open " << filename << " for writing" << std::endl;
    return;
  }
  
  // Write header
  csv << "Band,Corner,fMHz,";
  csv << "FundW,InsLossdB,FilterDissmW,";
  csv << "H2dBc,H3dBc,H4dBc,H5dBc,H6dBc,H7dBc,H8dBc,H9dBc,H10dBc,";
  csv << "C1IrmsA,C2IrmsA,C3IrmsA,C4IrmsA,L1IrmsA,L2IrmsA,L3IrmsA,";
  csv << "C1VpkV,C2VpkV,C3VpkV,C4VpkV,";
  csv << "C1Over,C2Over,C3Over,C4Over,L1Over,L2Over,L3Over" << std::endl;
  
  // Write data
  for (const auto& r : results) {
    csv << std::fixed << std::setprecision(6);
    csv << r.bandName << "," << r.cornerName << ",";
    csv << r.frequency / 1e6 << ",";
    csv << r.fundW << ",";
    csv << r.insLossdB << ",";
    csv << r.totalFiltDissW * 1000.0 << ",";
    
    // Harmonics 2-10
    for (int n = 2; n <= 10; ++n) {
      csv << r.harmdBc[n] << ",";
    }

    // Component stress
    csv << r.C1.Irms << "," << r.C2.Irms << "," << r.C3.Irms << "," << r.C4.Irms << ",";
    csv << r.L1.Irms << "," << r.L2.Irms << "," << r.L3.Irms << ",";
    csv << r.C1.Vpeak << "," << r.C2.Vpeak << "," << r.C3.Vpeak << "," << r.C4.Vpeak << ",";
    csv << (r.C1.overstressed ? "YES" : "NO") << ",";
    csv << (r.C2.overstressed ? "YES" : "NO") << ",";
    csv << (r.C3.overstressed ? "YES" : "NO") << ",";
    csv << (r.C4.overstressed ? "YES" : "NO") << ",";
    csv << (r.L1.overstressed ? "YES" : "NO") << ",";
    csv << (r.L2.overstressed ? "YES" : "NO") << ",";
    csv << (r.L3.overstressed ? "YES" : "NO");
    csv << std::endl;
  }
  
  csv.close();
  std::cout << "\nCSV exported to: " << filename << std::endl;
}

//=============================================================================
// Text Output Functions
//=============================================================================

void printResult(const HarmonicAnalyzer::BandResult& r) {
  std::cout << std::fixed << std::setprecision(3);

  std::cout << "  Frequency: " << r.frequency / 1e6 << " MHz" << std::endl;
  std::cout << "  Source Voltage (peak): " << r.vsourcePeak << " V" << std::endl;
  std::cout << "  DC Input Power: " << r.dcInputPower << " W" << std::endl;
  std::cout << "  RF to Filter: " << r.rfToFilter << " W" << std::endl;
  std::cout << "  Fundamental Power: " << r.fundW << " W" << std::endl;
  std::cout << "  Fundamental Filter Loss: " << r.fundDissipation * 1000.0 << " mW" << std::endl;
  std::cout << "  Insertion Loss: " << r.insLossdB << " dB" << std::endl;
  std::cout << "  Filter Dissipation: " << r.totalFiltDissW * 1000.0
            << " mW" << std::endl;
  
  std::cout << "\n  Harmonics (dBc):" << std::endl;
  for (int n = 2; n <= 10; ++n) {
    std::cout << "    H" << n << ": " << std::setw(8) << r.harmdBc[n] 
              << " dBc";
    if (n % 2 == 0) {
      std::cout << " (even - should be suppressed)";
    }
    std::cout << std::endl;
  }
  
  std::cout << "\n  Component Stress:" << std::endl;
  std::cout << "    C1: " << r.C1.Irms << " A(rms), "
            << r.C1.Vpeak << " V(pk)"
            << (r.C1.overstressed ? " **OVERSTRESSED**" : "") << std::endl;
  std::cout << "    C2: " << r.C2.Irms << " A(rms), "
            << r.C2.Vpeak << " V(pk)"
            << (r.C2.overstressed ? " **OVERSTRESSED**" : "") << std::endl;
  std::cout << "    C3: " << r.C3.Irms << " A(rms), "
            << r.C3.Vpeak << " V(pk)"
            << (r.C3.overstressed ? " **OVERSTRESSED**" : "") << std::endl;
  std::cout << "    C4: " << r.C4.Irms << " A(rms), "
            << r.C4.Vpeak << " V(pk)"
            << (r.C4.overstressed ? " **OVERSTRESSED**" : "") << std::endl;
  std::cout << "    L1: " << r.L1.Irms << " A(rms)"
            << (r.L1.overstressed ? " **OVERSTRESSED**" : "") << std::endl;
  std::cout << "    L2: " << r.L2.Irms << " A(rms)"
            << (r.L2.overstressed ? " **OVERSTRESSED**" : "") << std::endl;
  std::cout << "    L3: " << r.L3.Irms << " A(rms)"
            << (r.L3.overstressed ? " **OVERSTRESSED**" : "") << std::endl;
}

//=============================================================================
// Main Program
//=============================================================================

int main(int argc, char* argv[]) {
  // Parse command line argument for Veer voltage
  double veerVoltage = 55.5;  // Default value

  if (argc > 1) {
    try {
      veerVoltage = std::stod(argv[1]);
      if (veerVoltage < 0 || veerVoltage > 100) {
        std::cerr << "Error: Veer voltage must be between 0 and 100V" << std::endl;
        return 1;
      }
    } catch (const std::exception& e) {
      std::cerr << "Error: Invalid voltage argument. Usage: " << argv[0]
                << " [veer_voltage_in_volts]" << std::endl;
      return 1;
    }
  }

  std::cout << "========================================" << std::endl;
  std::cout << "NexRig LPF Harmonic Analysis" << std::endl;
  std::cout << "Veer Supply: " << veerVoltage << "V DC" << std::endl;
  std::cout << "========================================\n" << std::endl;
  
  // Band definitions - 7th-order Chebyshev LPF (0.25dB ripple, 200Ω impedance)
  // BandFilter(name, fLowMHz, fHighMHz, c1pF, c2pF, c3pF, c4pF, l1nH, l2nH, l3nH)
  // Prototype: g1=0.9309, g2=1.2243, g3=1.5827, g4=1.6896, g5=1.5827, g6=1.2243, g7=0.9309
  // C2 and C3 use 2 parallel caps for current sharing
  std::vector<BandFilter> bands = {
    // 160m: cutoff 2.4MHz
    BandFilter("160m", 1.8, 2.0,
               390, 330+360, 330+360, 390,
               18000, 20000, 18000),

    // 80m: cutoff 5.0MHz
    BandFilter("80m", 3.5, 4.0,
               180, 150+180, 150+180, 180,
               9100, 10000, 9100),

    // 60m: cutoff 6.5MHz
    BandFilter("60m", 5.3, 5.4,
               150, 160+110, 160+110, 150,
               6800, 7500, 6800),

    // 40m: cutoff 8.9MHz
    BandFilter("40m", 7.0, 7.3,
               110, 91+91, 91+91, 110,
               5100, 5600, 5100),

    // 30m: cutoff 11.5MHz
    BandFilter("30m", 10.1, 10.15,
               82, 75+75, 75+75, 82,
               3900, 4300, 3900),

    // 20m: cutoff 16MHz
    BandFilter("20m", 14.0, 14.35,
               56, 43+56, 43+56, 56,
               2700, 3000, 2700),

    // 17m/15m: cutoff 24MHz
    BandFilter("17m/15m", 18.068, 21.45,
               39, 33+36, 33+36, 39,
               1800, 2000, 1800),

    // 12m/10m: cutoff 32MHz
    BandFilter("12m/10m", 24.89, 29.7,
               30, 27+24, 27+24, 30,
               1500, 1600, 1500),
  };
  
  static const auto Ctol = 0.05;    // +/- 5%
  static const auto Ltol = 0.2;	    // +/- 20%

  // Tolerance corners
  std::vector<ToleranceCorner> corners = {
    {"Nominal",     1.00, 1.00, 1.00, 1.00},
    {"C+/L+",       1+Ctol, 1+Ltol, 0.90, 0.90},  // Caps high, inductors high
    {"C+/L-",       1+Ctol, 1-Ltol, 0.90, 1.10},  // Caps high, inductors low
    {"C-/L+",       1-Ctol, 1+Ltol, 1.10, 0.90},  // Caps low, inductors high
    {"C-/L-",       1-Ctol, 1-Ltol, 1.10, 1.10},  // Caps low, inductors low
    {"WorstLoss",   1.00, 1.00, 1.30, 1.30}   // Nominal values, high loss
  };
  
  // Baseline parasitics (nominal)
  const double ESRbase = 0.125;  // Ohms (C0G at HF)
  const double DCRbase = 0.055;  // Ohms (Coilcraft 132-xx)
  
  // Square wave source - Veer supply voltage from command line (default 55.5V)
  // Transformer: 3:7 turns (36:196Ω impedance ratio)
  SquareWaveSource source(veerVoltage, 3e-9, 3e-9);
  
  // Storage for all results
  std::vector<HarmonicAnalyzer::BandResult> allResults;
  
  // Analyze each band
  for (const auto& band : bands) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Band: " << band.name << std::endl;
    std::cout << "========================================" << std::endl;

    // Analyze each corner
    for (const auto& corner : corners) {
      std::cout << "\n--- Corner: " << corner.name << " ---" << std::endl;

      // Apply tolerances
      double C1 = band.c1 * corner.Cfactor;
      double C2 = band.c2 * corner.Cfactor;
      double C3 = band.c3 * corner.Cfactor;
      double C4 = band.c4 * corner.Cfactor;
      double L1 = band.l1 * corner.Lfactor;
      double L2 = band.l2 * corner.Lfactor;
      double L3 = band.l3 * corner.Lfactor;
      double ESR = ESRbase * corner.ESRfactor;
      double DCR = DCRbase * corner.DCRfactor;

      // Create filter with corner values (7th-order)
      ChebyshevLPF filter(C1, C2, C3, C4, L1, L2, L3, ESR, DCR);

      // For combined bands, analyze each sub-band separately
      if (band.name == "17m/15m") {
        // 17m band: 18.068-18.168 MHz
        std::cout << "\n[17m Bottom]" << std::endl;
        auto m17low = HarmonicAnalyzer::analyze(
          filter, source, "17m", corner.name, 18.068e6, 10
        );
        printResult(m17low);
        allResults.push_back(m17low);

        std::cout << "\n[17m Top]" << std::endl;
        auto m17high = HarmonicAnalyzer::analyze(
          filter, source, "17m", corner.name, 18.168e6, 10
        );
        printResult(m17high);
        allResults.push_back(m17high);

        // 15m band: 21.0-21.45 MHz
        std::cout << "\n[15m Bottom]" << std::endl;
        auto m15low = HarmonicAnalyzer::analyze(
          filter, source, "15m", corner.name, 21.0e6, 10
        );
        printResult(m15low);
        allResults.push_back(m15low);

        std::cout << "\n[15m Top]" << std::endl;
        auto m15high = HarmonicAnalyzer::analyze(
          filter, source, "15m", corner.name, 21.45e6, 10
        );
        printResult(m15high);
        allResults.push_back(m15high);
      }
      else if (band.name == "12m/10m") {
        // 12m band: 24.89-24.99 MHz
        std::cout << "\n[12m Bottom]" << std::endl;
        auto m12low = HarmonicAnalyzer::analyze(
          filter, source, "12m", corner.name, 24.89e6, 10
        );
        printResult(m12low);
        allResults.push_back(m12low);

        std::cout << "\n[12m Top]" << std::endl;
        auto m12high = HarmonicAnalyzer::analyze(
          filter, source, "12m", corner.name, 24.99e6, 10
        );
        printResult(m12high);
        allResults.push_back(m12high);

        // 10m band: 28.0-29.7 MHz
        std::cout << "\n[10m Bottom]" << std::endl;
        auto m10low = HarmonicAnalyzer::analyze(
          filter, source, "10m", corner.name, 28.0e6, 10
        );
        printResult(m10low);
        allResults.push_back(m10low);

        std::cout << "\n[10m Top]" << std::endl;
        auto m10high = HarmonicAnalyzer::analyze(
          filter, source, "10m", corner.name, 29.7e6, 10
        );
        printResult(m10high);
        allResults.push_back(m10high);
      }
      else {
        // Regular bands - analyze bottom and top
        std::cout << "\n[Bottom of Band]" << std::endl;
        auto low = HarmonicAnalyzer::analyze(
          filter, source, band.name, corner.name, band.fLow, 10
        );
        printResult(low);
        allResults.push_back(low);

        std::cout << "\n[Top of Band]" << std::endl;
        auto high = HarmonicAnalyzer::analyze(
          filter, source, band.name, corner.name, band.fHigh, 10
        );
        printResult(high);
        allResults.push_back(high);
      }
    }
  }
  
  // Export to CSV
  exportResultsToCSV(allResults, "harmonics.csv");
  
  std::cout << "\n========================================" << std::endl;
  std::cout << "Analysis Complete" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "\nTotal analyses performed: " << allResults.size() << std::endl;
  std::cout << "Results saved to: harmonics.csv" << std::endl;
  
  return 0;
}
