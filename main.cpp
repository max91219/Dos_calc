#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/assert.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

namespace mpi = boost::mpi;

unsigned int good_seed(mpi::communicator c);
double disp(double kx, double ky, double kz, bool use_tp);
int find_bin(double val, double bottom, double bin_width);

int main(int argc, char* argv[]) {

  // Set up the program options
  namespace po = boost::program_options;
  po::options_description description("Options");
  description.add_options()
    ("o_file,o", po::value<std::string>()->default_value("dos.dat"), "output file name")
    ("use_tp,t", po::bool_switch()->default_value(false), "use tp equals 1/4")
    ("n_sample,n", po::value<long int>()->default_value(2000000), "number of samples per thread")
    ("num_vals,N", po::value<int>()->default_value(200000), "number of points in the epsilon grid");

  po::variables_map vm;

  try {
    po::store(po::command_line_parser(argc, argv).options(description).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << "USAGE: " << argv[0] << std::endl;
      std::cout << description << std::endl;
      return 0;
    }
  } catch (po::error &e) {
    std::cerr << e.what() << std::endl;
    std::cout << "USAGE: " << argv[0] << std::endl;
    std::cout << description << std::endl;
    return EXIT_FAILURE;
  }

  int N = vm["num_vals"].as<int>();
  long int n = vm["n_sample"].as<long int>();
  std::string o_file = vm["o_file"].as<std::string>();
  bool use_tp = vm["use_tp"].as<bool>();

  // Set up MPI
  mpi::environment env(argc, argv);
  mpi::communicator comm;
  int rank = comm.rank();

  //Set up the RNG
  boost::random::mt19937 RNG_engine;

  double upper_band_edge = use_tp ? (4.0 - 2.0*(0.25)) : 4.0;
  double lower_band_edge = use_tp ? upper_band_edge - (16.0 + 1.0)  : upper_band_edge - 16.0;

  if (rank == 0) std::cout << "lower band edge: " << lower_band_edge << std::endl;
  if (rank == 0) std::cout << "upper band edge: " << upper_band_edge << std::endl;

  std::vector<double> eps_vals(N+1);
  std::vector<double> dos(N+1);
  double d_eps = (upper_band_edge - lower_band_edge) / (double) N;
  double d_eps_2 = d_eps/2.0;
  double kx, ky, kz, disp_val;

  if (rank == 0) std::cout << "N: " << N << std::endl;
  if (rank == 0) std::cout << "n: " << n << std::endl;
  if (rank == 0) std::cout << "d_eps: " << d_eps << std::endl;
  if (rank == 0) std::cout << "d_eps_2: " << d_eps_2 << std::endl;
  if (rank == 0) std::cout << "o_file: " << o_file << std::endl;
  if (rank == 0) std::cout << "use_tp: " << use_tp << std::endl;

  for (int i = 0; i<eps_vals.size() ; i++) {
    eps_vals[i]= lower_band_edge + (i * d_eps);
    dos[i] = 0.0;
    //std::cout << eps_vals[i] << std::endl;
  }

  //Seed the RNG
  RNG_engine.seed(good_seed(comm));
  boost::uniform_real<> uni_dist_0_2PI(0,2.0*boost::math::constants::pi<double>());
  boost::variate_generator<boost::random::mt19937&, boost::uniform_real<> > rng(RNG_engine, uni_dist_0_2PI);

  for (long int j = 0; j <= n; j++) {

    if(j%(n/10) == 0 && rank == 0) std::cout << "on sample: " << j << std::endl;

    // Generate random wavevector
    kx = rng();
    BOOST_ASSERT_MSG(0.0 <= kx <= 2.0 * boost::math::constants::pi<double>(), "kx out of range [0,2Pi]");
    ky = rng();
    BOOST_ASSERT_MSG(0.0 <= ky <= 2.0 * boost::math::constants::pi<double>(), "ky out of range [0,2Pi]");
    kz = rng();
    BOOST_ASSERT_MSG(0.0 <= kz <= 2.0 * boost::math::constants::pi<double>(), "kz out of range [0,2Pi]");

    disp_val = disp(kx, ky, kz, use_tp);
    BOOST_ASSERT_MSG(disp_val >= lower_band_edge, "disp below lower band edge");
    BOOST_ASSERT_MSG(disp_val <= upper_band_edge, "disp above upper band edge");

    int found_i = find_bin(disp_val, lower_band_edge, d_eps);
    BOOST_ASSERT_MSG(0<= found_i <= eps_vals.size(), "Mapped to incorrect index");

    dos[found_i]+= 1.0/d_eps;
    //std::cout << "rank: " << rank << " val : " << disp_val << " eps: " <<  eps_vals[found_i] << std::endl;
    //std::cout << "rank : " << rank << " disp : " << disp_val << std::endl;
  }

  if (rank == 0) std::cout << "Done Calculating." << std::endl;

  if (rank == 0) {
    std::vector<double> dos_res(N+1);
    //Reduce the data
    mpi::reduce(comm, dos, dos_res, std::plus<double>(), 0);
    double n_procs = comm.size();
    double val = 0.0;

    std::ofstream ofile;
    ofile.open(o_file.c_str());
    if (ofile.is_open()){
      for (int i = 0; i < dos_res.size(); i++){
        val = dos_res[i];
        val /= n_procs;
        val /= n;

        ofile << eps_vals[i] << " " << val << std::endl;
      }
      ofile.close();
    }
  } else {
    //Reduce the data
    mpi::reduce(comm, dos, std::plus<double>(), 0);

  }


  return 0;
}

int find_bin(double val, double bottom, double bin_width) {
  int i = static_cast<int>((val - bottom)/bin_width);
  return fmod((val - bottom), bin_width) >= 0.5 ? i + 1 : i;
}

unsigned int good_seed(mpi::communicator c) {
  unsigned int random_seed, random_seed_a, random_seed_b;

  std::ifstream file ("/dev/urandom", std::ios::binary);
  if (file.is_open()) {
    char * memblock;
    int size = sizeof(int);
    memblock = new char [size];
    file.read (memblock, size);
    file.close();
    random_seed_a = *reinterpret_cast<int*>(memblock);
    delete[] memblock;
  }

  random_seed_b = std::time(0);
  random_seed = random_seed_a xor random_seed_b xor c.rank();
  std::cout << "random_seed =  " << random_seed << "  rank = " << c.rank() << std::endl;
  return random_seed;
}

double disp(double kx, double ky, double kz, bool use_tp) {
  double dispersion = 0.0;

  dispersion = -4.0 * (std::cos(kx/2.0)*std::cos(ky/2.0) + std::cos(kz/2.0)*std::cos(ky/2.0) + std::cos(kx/2.0)*std::cos(kz/2.0));

  if (use_tp) {
    dispersion -= 2.0 * 0.25 * (std::cos(kx) + std::cos (ky) + std::cos(kz));
  }

  return dispersion;
}
