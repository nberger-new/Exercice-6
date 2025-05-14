#include "ConfigFile.tpp"
#include <chrono>
#include <cmath>
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
    if (x_ >= xL_ && x_ <= xa_) { V_ = 0.5*pow(omega_0, 2)*pow((x_-xa_)/(xL_ - xa_), 2); }
    else if (x_ >= xL_ && x_ <= xa_) { V_ = 0.5*pow(std::sin(M_PI*(x_ - xa_)/(xb_ - xa_)), 2) }
    else if (x_ >= xL_ && x_ <= xa_) { V_ = 0.5*pow(omega_0, 2)*pow((x_-xb_)/(xR_ - xb_), 2); }
    else { cerr << "Choisir une position valide !" << endl;}
    return V_;
}

// Declaration des diagnostiques de la particule d'apres sa fonction d'onde psi :
//  - prob: calcule la probabilite de trouver la particule dans un intervalle [x_i, x_j]
//  - E:    calcule son energie,
//  - xmoy: calcule sa position moyenne,
//  - x2moy:calcule sa position au carre moyenne,
//  - pmoy: calcule sa quantite de mouvement moyenne,
//  - p2moy:calcule sa quantite de mouvement au carre moyenne.


// TODO: calculer la probabilite de trouver la particule dans un intervalle [x_i, x_j]
double prob(int I, int J, vec_cmplx Psi)
{
    sum  = 0;
    for (i = I; i < J; i++) {
        sum += 0.5*(Psi[i]*std::conj(Psi[i]) + Psi[i+1]*std::conj(Psi[i+1]))
    }
    return sum;
}

// TODO calculer l'energie
double E(vec_cmplx Psi){

    return 0;
}

// TODO calculer xmoyenne
double xmoy()
{
    return 0;
}

// TODO calculer x.^2 moyenne
double x2moy()
{
    return 0;
}

// TODO calculer p moyenne
double pmoy()
{
    return 0;
}

// TODO calculer p.^2 moyenne
double p2moy(double h, vector <double> Psi, vector <double> upper_, vector <double> diag_)
{
    double sum = 0;
    int N = Psi.size();
    for (int i = 0; i < N - 1; ++i) {
        double p_term_i = -1/h*(Psi[i -1] - 2*Psi[i] + Psi[i+1]);
        double p_term_i_next = -1/h*(Psi[i] - 2*Psi[i + 1] + Psi[i+2]);
        sum = Psi[i]*p_term_i + Psi[i+1]*p_term_i_next;
    }
    return 0.5*sum;
}

// TODO calculer la normalization
vec_cmplx normalize(vec_cmplx const& Psi, double const& dx)
{
    vec_cmplx psi_norm(psi.size(), 0.);
    int N = Psi.size();
    double sum = 0.;
    for (int i = 0; i < N - 1; ++i) {
        sum += std::conj(Psi[i])*Psi[i] + std::conj(Psi[i+1])*Psi[i+1];
    }
    sum *= dx/2;
    for (int i = 0; i < N; ++i) {
        psi_norm[i] = Psi[i]/std::sqrt(sum);
    }
    return psi_norm;
}




int
main(int argc, char** argv)
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
    double k0 = 2*PI*n/L;

    int Npoints = Nintervals + 1;
    double dx = (xR - xL) / Nintervals;
    double dt = tfin / Nsteps;

    const auto simulationStart = std::chrono::steady_clock::now();

    // Maillage :
    vector<double> x(Npoints);
    for (int i(0); i < Npoints; ++i)
        x[i] = xL + i * dx;

    // Initialisation de la fonction d'onde :
    vec_cmplx psi(Npoints, 0.);

    // initialization time and position to check Probability
    double t = 0;
    unsigned int Nx0 = floor((0 - xL)/(xR-xL)*Npoints); //chosen xR*0.5 since top of potential is at half x domain
  
    // TODO initialize psi
    for (int i(0); i < Npoints; ++i)
    	psi[i] = std::exp(complex_i*k0*x[i])*std::exp(-pow((x[i]-x[0]), 2)/(2*sigma0*sigma0));
   
    // Modifications des valeurs aux bords :
    psi[0] = complex<double>(0., 0.);
    psi[Npoints - 1] = complex<double>(0., 0.);
    
    // Normalisation :
    psi = normalize(psi, dx);

    // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
    vec_cmplx dH(Npoints), aH(Nintervals), cH(Nintervals); // matrice Hamiltonienne
    vec_cmplx dA(Npoints), aA(Nintervals),
      cA(Nintervals); // matrice du membre de gauche de l'equation (4.100)
    vec_cmplx dB(Npoints), aB(Nintervals),
      cB(Nintervals); // matrice du membre de droite de l'equation (4.100)

    complex<double> a = complex_i * hbar * dt / (4.*m*dx*dx); // Coefficient complexe a de l'equation (4.100)
    complex<double> b = complex_i * dt / (2*hbar); // Coefficient complexe a de l'equation (4.100)

    // TODO: calculer les éléments des matrices A, B et H.
    // Ces matrices sont stockées sous forme tridiagonale, d:diagonale, c et a: diagonales
    // supérieures et inférieures
    for (int i(0); i < Npoints; ++i) // Boucle sur les points de maillage
    {
        V_i = V(V0, x[i], xL, xR, xa, xb, om0);
        dH[i] = -hbar*hbar/(4*dx*dx) + V_i;
        dA[i] = 1 + 2*a + b*V_i;
        dB[i] = 1 + 2*a + b*V_i;
    }
    for (int i(0); i < Nintervals; ++i) // Boucle sur les intervalles
    {
        aH[i] = -hbar*hbar/(4*dx*dx);
        aA[i] = -a;
        aB[i] = a;
        cH[i] = -hbar*hbar/(4*dx*dx);
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

    // Fichiers de sortie :
    string output = configFile.get<string>("output");

    ofstream fichier_potentiel((output + "_pot.out").c_str());
    fichier_potentiel.precision(15);
    for (int i(0); i < Npoints; ++i)
        fichier_potentiel << x[i] << " " <<V() << endl;
    fichier_potentiel.close();

    ofstream fichier_psi((output + "_psi2.out").c_str());
    fichier_psi.precision(6);

    ofstream fichier_observables((output + "_obs.out").c_str());
    fichier_observables.precision(15);

    // t0 writing
    for (int i(0); i < Npoints; ++i){
        fichier_psi << abs(psi[i])  << " " << real(psi[i]) << " "  << imag(psi[i]) << " ";
        }
    fichier_psi << endl;

    // Ecriture des observables :
    // TODO: introduire les arguments des fonctions prob, E, xmoy, x2moy, pmoy et p2moy
    //       en accord avec la façon dont vous les aurez programmés plus haut
    fichier_observables << t << " " << prob() << " " << prob()
                << " " << E() << " " << xmoy () << " "  
                << x2moy() << " " << pmoy () << " " << p2moy() << endl; 

    // Boucle temporelle :    
    while (t < tfin) {

        // Multiplication psi_tmp = B * psi :
        vec_cmplx psi_tmp(Npoints, 0.);
        for (int i(0); i < Npoints; ++i)
            psi_tmp[i] = dB[i] * psi[i];
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
        fichier_observables << t << " " << prob() << " " << prob()
                    << " " << E() << " " << xmoy () << " "  
                    << x2moy() << " " << pmoy () << " " << p2moy() << endl; 

    } // Fin de la boucle temporelle





    fichier_observables.close();
    fichier_psi.close();

    const auto simulationEnd = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsedSeconds = simulationEnd - simulationStart;
    std::cout << "Simulation finished in " << setprecision(3) << elapsedSeconds.count()
              << " seconds" << std::endl;
}
