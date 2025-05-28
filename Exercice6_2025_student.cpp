#include "ConfigFile.tpp"
#include <chrono>
#include <cmath>
#include <algorithm>
#include <complex> // Pour les nombres complexes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
typedef vector<complex<double>> vec_cmplx;

// Fonction resolvant le systeme d'equations A * solution = rhs
// où A est une matrice tridiagonale
template<class T>
void triangular_solve(vector<T> const& diag,  vector<T> const& lower, vector<T> const& upper,
                 vector<T> const& rhs, vector<T>& solution)
{
    vector<T> new_diag = diag;
    vector<T> new_rhs = rhs;

    // forward elimination
    for (int i(1); i < diag.size(); ++i) {
        T pivot = lower[i - 1] / new_diag[i - 1];
        new_diag[i] -= pivot * upper[i - 1];
        new_rhs[i] -= pivot * new_rhs[i - 1];
    }

    solution.resize(diag.size());

    // solve last equation
    solution[diag.size() - 1] = new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

    // backward substitution
    for (int i = diag.size() - 2; i >= 0; --i) {
        solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];
    }
}

// TODO Potentiel V(x) :
double V(double V0_, double x_, double xL_, double xR_, double xa_, double xb_, double omega_0)
{
    double V_(0);
    if (x_ >= xL_ && x_ < xa_) { V_ = 0.5*pow(omega_0, 2)*pow((x_-xa_)/(1 - xa_/xL_), 2); }
    else if (x_ >=  xa_ && x_ <= xb_) {
        if (xb_ - xa_ != 0){ V_ = V0_*pow(std::sin(M_PI*(x_ - xa_)/(xb_ - xa_)), 2); }
        else { V_ = V0_*pow(std::sin(M_PI*(x_ - xa_)), 2); }
        //cout << "V_ = " << V_ << endl;

    }
    else if (x_ > xb_ && x_ <= xR_) { V_ = 0.5*pow(omega_0, 2)*pow((x_-xb_)/(1 - xb_/xR_), 2); }
    else { cerr << "Choisir une position valide !" << endl;}
    return V_;
}

// if (xb_ - xa_ != 0){ V_ = V0*pow(std::sin(M_PI*(x_ - xa_)/(xb_ - xa_)), 2);}
// else {V_ = V0_*pow(std::sin(M_PI*(x_ - xa_)), 2); }

// Declaration des diagnostiques de la particule d'apres sa fonction d'onde psi :
//  - prob: calcule la probabilite de trouver la particule dans un intervalle [x_i, x_j]
//  - E:    calcule son energie,
//  - xmoy: calcule sa position moyenne,
//  - x2moy:calcule sa position au carre moyenne,
//  - pmoy: calcule sa quantite de mouvement moyenne,
//  - p2moy:calcule sa quantite de mouvement au carre moyenne.


// TODO: calculer la probabilite de trouver la particule dans un intervalle [x_i, x_j]
double prob(int I, int J, vec_cmplx const& Psi, double const& dx)
{
    double sum  = 0;
    for (int i = I; i < J; i++) {
        //cout << std::norm(Psi[i]*conj(Psi[i])) << endl;
        sum += 0.5*dx*(std::norm(Psi[i]) + std::norm(Psi[i+1]));
    }
    return sum;
}

// TODO calculer l'energie
double E(vec_cmplx const& Psi, vec_cmplx const& H_up, vec_cmplx const& H_down, vec_cmplx const& H_diag, double const& dx)
{
    double E = 0;
    vec_cmplx HPsi(Psi.size(), 0.0);
    for (int i = 1; i < Psi.size()-1; ++i) {
        HPsi[i] += H_diag[i]*Psi[i];
        HPsi[i] += H_up[i]*Psi[i+1];
        HPsi[i+1] += H_down[i]*Psi[i];
        //if (i > 0) { HPsi[i] += H_down[i-1]*Psi[i-1]; }
        //HPsi[i] += H_up[i+1]*Psi[i+1];
    }
    HPsi[Psi.size()-1] = 0.0;
    for (int i = 0; i < Psi.size()-1; ++i) {
        E += std::real(std::conj(Psi[i])*HPsi[i]) + std::real(std::conj(Psi[i+1])*HPsi[i+1]);
    }
    E *= dx/(2.);
    return E;
}

// E += std::real(std::conj(Psi[i])*H_diag[i]*Psi[i]);
// if (i > 0) { E += std::real(std::conj(Psi[i])*H_up[i-1]*Psi[i]); }
// if (i < Psi.size() - 1) { E += std::real(std::conj(Psi[i])*H_down[i+1]*Psi[i]  + std::conj(Psi[i+1])*H_diag[i+1]*Psi[i+1] + std::conj(Psi[i+1])*H_up[i]*Psi[i+1]); }
// if (i < Psi.size() - 2) { E += std::real(std::conj(Psi[i+1])*H_down[i+2]*Psi[i+1]); }
// cout << " At i = " << i << " E = " << E << endl;

// TODO calculer xmoyenne
double xmoy(vec_cmplx const& Psi, vector<double> const& x_, double const& dx)
{
    double x_mean = 0;
    for (int i = 0; i < Psi.size()-1; ++i) {
        x_mean += std::real(std::conj(Psi[i])*x_[i]*Psi[i] + std::conj(Psi[i+1])*x_[i+1]*Psi[i+1]);
    }
    return dx*x_mean/2.;
}

// TODO calculer x.^2 moyenne
double x2moy(vec_cmplx const& Psi, vector<double> const& x_, double const& dx)
{
    double x2_mean = 0;
    for (int i = 0; i < Psi.size(); ++i) {
        x2_mean += std::real(std::conj(Psi[i])*x_[i]*x_[i]*Psi[i] + std::conj(Psi[i+1])*x_[i+1]*x_[i+1]*Psi[i+1]);
    }
    return dx*x2_mean/2.;
}

// TODO calculer p moyenne
double pmoy(vec_cmplx const& Psi, double const& dx)
{
    double p_mean = 0;
    complex<double> complex_i = complex<double>(0, 1);

    vec_cmplx der_Psi(Psi.size(), 0.);

    der_Psi[0] = Psi[1]/dx;
    for (int i = 1; i < Psi.size() - 2; ++i) {
        der_Psi[i] = (Psi[i+1] - Psi[i - 1])/(2.*dx);
    }
    der_Psi[Psi.size() - 1] = Psi[Psi.size()-2]/dx;

    for (int i = 1; i < Psi.size() - 1; ++i) {
        p_mean += std::real(-complex_i*(std::conj(Psi[i])*der_Psi[i] + std::conj(Psi[i+1])*der_Psi[i+1]));
    }
    //cout << "p_mean = " << p_mean << endl;
    return std::real(dx*p_mean/2.);
}

// TODO calculer p.^2 moyenne
double p2moy(vec_cmplx const& Psi, double const& dx)
{
    double sum = 0;
    vec_cmplx der_2Psi(Psi.size(), 0.);
    for (int i = 1; i < Psi.size() - 2; ++i) {
        der_2Psi[i] = (Psi[i+1] - Psi[i] - Psi[i] + Psi[i-1])/(dx*dx);
    }
    for (int i = 0; i < Psi.size() - 1; ++i) {
        sum += std::real(std::conj(Psi[i])*der_2Psi[i] + std::conj(Psi[i+1])*der_2Psi[i+1]);
    }
    return -dx*sum/2.;
}

// TODO calculer la normalization
vec_cmplx normalize(vec_cmplx const& Psi, double const& dx)
{
    vec_cmplx psi_norm(Psi.size(), 0.);
    double sum = 0.;
    for (int i = 0; i < Psi.size() - 1; ++i) {
        //cout << std::real(std::conj(Psi[i])*Psi[i] + std::conj(Psi[i+1])*Psi[i+1]) << " =? " << std::real(norm(Psi[i]) + norm(Psi[i+1])) << endl;
        sum += 0.5*dx*std::real(std::conj(Psi[i])*Psi[i] + std::conj(Psi[i+1])*Psi[i+1]);
    }
    for (int i = 0; i < Psi.size(); ++i) {
        psi_norm[i] = Psi[i]/std::sqrt(sum);
    }
    // sum = 0.;
    // for (int i = 0; i < Psi.size() - 1; ++i) {
    //     //cout << std::real(std::conj(Psi[i])*Psi[i] + std::conj(Psi[i+1])*Psi[i+1]) << " =? " << std::real(norm(Psi[i]) + norm(Psi[i+1])) << endl;
    //     sum += 0.5*dx*std::real(std::norm(psi_norm[i]) + std::norm(psi_norm[i+1]));
    //     cout << "sum = " << sum << endl;
    // }
    return psi_norm;
}




int main(int argc, char** argv)
{
    complex<double> complex_i = complex<double>(0, 1); // Nombre imaginaire i
    const double PI = 3.1415926535897932384626433832795028841971e0;

    string inputPath("configuration.in.example"); // Fichier d'input par defaut
    if (argc > 1) // Fichier d'input specifie par l'utilisateur ("./Exercice8 config_perso.in")
        inputPath = argv[1];

    ConfigFile configFile(
      inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for (int i(2); i < argc;
         ++i) // Input complementaires ("./Exercice8 config_perso.in input_scan=[valeur]")
        configFile.process(argv[i]);

    // Set verbosity level. Set to 0 to reduce printouts in console.
    const int verbose = configFile.get<int>("verbose");
    configFile.setVerbosity(verbose);

    // Parametres physiques :
    double hbar = 1.;
    double m = 1.;
    double tfin = configFile.get<double>("tfin");
    double xL = configFile.get<double>("xL");
    double xR = configFile.get<double>("xR");
    double xa = configFile.get<double>("xa");
    double xb = configFile.get<double>("xb");
    double V0 = configFile.get<double>("V0");
    double om0 = configFile.get<double>("om0");
    double n  = configFile.get<int>("n"); // Read mode number as integer, convert to double

    double x0 = configFile.get<double>("x0");
    double sigma0 = configFile.get<double>("sigma_norm") * (xR - xL);

    int Nsteps = configFile.get<int>("Nsteps");
    int Nintervals = configFile.get<int>("Nintervals");

    // TODO: initialiser le paquet d'onde, equation (4.116) du cours
    double k0 = 2*PI*n/(xR - xL);

    int Npoints = Nintervals + 1;
    double dx = (xR - xL) / Nintervals;
    double dt = tfin / Nsteps;

    const auto simulationStart = std::chrono::steady_clock::now();

    // Maillage :
    vector<double> x(Npoints);
    for (int i(0); i < Npoints; ++i) {
        x[i] = xL + i * dx;
        //cout << " i = " << i << " x[i] = " << x[i] << endl;
    }
    // Initialisation de la fonction d'onde :
    vec_cmplx psi(Npoints, 0.);

    // initialization time and position to check Probability
    double t = 0;
    unsigned int Nx0 = floor((0 - xL)/(xR-xL)*Npoints); //chosen xR*0.5 since top of potential is at half x domain
  
    // TODO initialize psi
    for (int i(0); i < Npoints; ++i) {
    	psi[i] = std::exp(complex_i*k0*x[i])*std::exp(-pow((x[i]-x0), 2)/(2*sigma0*sigma0));
      // cout << "i = " << i << " with " << psi[i] << endl;
    }
   
    // Modifications des valeurs aux bords :
    psi[0] = complex<double>(0., 0.);
    psi[Npoints - 1] = complex<double>(0., 0.);
    
    // Normalisation :
    psi = normalize(psi, dx);
    for (int i(0); i < Npoints; ++i) {
      //cout << "i = " << i << " with " << psi[i] << endl;
    }
    // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
    vec_cmplx dH(Npoints), aH(Nintervals), cH(Nintervals); // matrice Hamiltonienne
    vec_cmplx dA(Npoints), aA(Nintervals),
      cA(Nintervals); // matrice du membre de gauche de l'equation (4.100)
    vec_cmplx dB(Npoints), aB(Nintervals),
      cB(Nintervals); // matrice du membre de droite de l'equation (4.100)

    complex<double> a = complex_i * hbar * dt / (4.*m*dx*dx); // Coefficient complexe a de l'equation (4.100)
    //cout << " a = " << a << endl;
    complex<double> b = complex_i * dt / (2.*hbar); // Coefficient complexe a de l'equation (4.100)
    //cout << " b = " << b << endl;
    complex<double> complex_1 = complex<double>(1., 0);
    // TODO: calculer les éléments des matrices A, B et H.
    // Ces matrices sont stockées sous forme tridiagonale, d:diagonale, c et a: diagonales
    // supérieures et inférieures
    for (int i(0); i < Npoints; ++i) // Boucle sur les points de maillage
    {
        double V_i = V(V0, x[i], xL, xR, xa, xb, om0);
        //cout << " V[i] = " << V_i << endl;
        dH[i] = hbar*hbar/(m*dx*dx) + V_i;
        dA[i] = complex_1 + 2.*a + b*V_i;
        dB[i] = complex_1 - 2.*a - b*V_i;
    }
    for (int i(0); i < Nintervals; ++i) // Boucle sur les intervalles
    {
        aH[i] = -hbar*hbar/(2.*m*dx*dx);
        aA[i] = -a;
        aB[i] = a;
        cH[i] = -hbar*hbar/(2.*m*dx*dx);
        cA[i] = -a;
        cB[i] = a;
    }

    // Conditions aux limites: psi nulle aux deux bords
    // TODO: Modifier les matrices A et B pour satisfaire les conditions aux limites
    aA[0] = 0;
    aB[0] = 0;
    cA[0] = 0;
    cB[0] = 0;
    aA[Nintervals-1] = 0;
    aB[Nintervals-1] = 0;
    cA[Nintervals-1] = 0;
    cB[Nintervals-1] = 0;


    // for (int i(0); i < Npoints; ++i) { cout << " i = " << i << " dH[i] = " << dH[i] << " dA[i] = " << dA[i] << " dB[i] = " << dB[i] << endl; }
    // for (int i(0); i < Nintervals; ++i) {
    //     cout << " i = " << i << " aH[i] = " << aH[i] << " aA[i] = " << aA[i]  << " aB[i] = " << aB[i] << " cH[i] = " << cH[i] << " cA[i] = " << cA[i] << " cB[i] = " << cB[i] << endl;
    // }

    // Fichiers de sortie :
    string output = configFile.get<string>("output");

    ofstream fichier_potentiel((output + "_pot.out").c_str());
    fichier_potentiel.precision(15);
    for (int i(0); i < Npoints; ++i)
        fichier_potentiel << x[i] << " " << V(V0, x[i], xL, xR, xa, xb, om0) << endl;
    fichier_potentiel.close();

    ofstream fichier_psi((output + "_psi.out").c_str());
    fichier_psi.precision(6);

    ofstream fichier_observables((output + "_obs.out").c_str());
    fichier_observables.precision(15);

    ofstream fichier_inc((output + "_inc.out").c_str());
    fichier_inc.precision(15);

    // t0 writing
    for (int i(0); i < Npoints; ++i){
        fichier_psi << abs(psi[i])  << " " << real(psi[i]) << " "  << imag(psi[i]) << " ";
        }
    fichier_psi << endl;

    // Ecriture des observables :
    // TODO: introduire les arguments des fonctions prob, E, xmoy, x2moy, pmoy et p2moy
    //       en accord avec la façon dont vous les aurez programmés plus haut
    double target = 0;
    double epsilon = 0.01;
    auto it = find_if(x.begin(), x.end(), [=](double b) {
        return std::abs(target - b) < epsilon;
    });
    if (it == x.end()) { cerr << "Lower tolerance !!" << endl;}
    else {
    fichier_observables << t << " " << prob(0, it - x.begin(), psi, dx) << " " << prob(it - x.begin(), psi.size()-1, psi, dx)
                << " " << E(psi, cH, aH, dH, dx) << " " << xmoy (psi, x, dx) << " "
                << x2moy(psi, x, dx) << " " << pmoy (psi, dx) << " " << p2moy(psi, dx) << endl;
    }

    double delta_mean_x = std::sqrt(x2moy(psi, x, dx) - pow(xmoy (psi, x, dx), 2));
    double delta_mean_p = std::sqrt(p2moy(psi, dx) - pow(pmoy (psi, dx), 2));

    fichier_inc << t << " " << delta_mean_x << " " << delta_mean_p << endl;

    // Boucle temporelle :    
    while (t < tfin) {

        // Multiplication psi_tmp = B * psi :
        vec_cmplx psi_tmp(Npoints, 0.);
        for (int i(0); i < Npoints; ++i) {
            psi_tmp[i] = dB[i] * psi[i];
        }
        for (int i(0); i < Nintervals; ++i) {
            psi_tmp[i] += cB[i] * psi[i + 1];
            psi_tmp[i + 1] += aB[i] * psi[i];
        }

        // Resolution de A * psi = psi_tmp :
        triangular_solve(dA, aA, cA, psi_tmp, psi);
        t += dt;

        // t0 writing
        for (int i(0); i < Npoints; ++i){
            fichier_psi << abs(psi[i])  << " " << real(psi[i]) << " "  << imag(psi[i]) << " ";
            }
        fichier_psi << endl;

        // Ecriture des observables :
	// TODO: introduire les arguments des fonctions prob, E, xmoy, x2moy, pmoy et p2moy
	//       en accord avec la façon dont vous les aurez programmés plus haut
        fichier_observables << t << " " << prob(0, it - x.begin(), psi, dx) << " " << prob(it - x.begin(), psi.size()-1, psi, dx)
                            << " " << E(psi, cH, aH, dH, dx) << " " << xmoy (psi, x, dx) << " "
                            << x2moy(psi, x, dx) << " " << pmoy (psi, dx) << " " << p2moy(psi, dx) << endl;

        double delta_mean_x = std::sqrt(x2moy(psi, x, dx) - pow(xmoy (psi, x, dx), 2));
        double delta_mean_p = std::sqrt(p2moy(psi, dx) - pow(pmoy (psi, dx), 2));

        fichier_inc << t << " " << delta_mean_x << " " << delta_mean_p << endl;

    } // Fin de la boucle temporelle

    fichier_observables.close();
    fichier_psi.close();
    fichier_inc.close();

    const auto simulationEnd = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsedSeconds = simulationEnd - simulationStart;
    std::cout << "Simulation finished in " << setprecision(3) << elapsedSeconds.count()
              << " seconds" << std::endl;
}
