#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <boost/program_options.hpp>

int main(int ac, char** av) {

  namespace po = boost::program_options;
  po::options_description description("Options");
  description.add_options()
    ("out_file,o", po::value<std::string>()->required(), "Name of output HDF5 file")
    ("in_files,i", po::value<std::vector<std::string> >()->multitoken()->required(), "input files");
  po::variables_map vm;

  try {
    po::store(po::command_line_parser(ac, av).options(description).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << "USAGE: " << av[0] << std::endl;
      std::cout << description << std::endl;
      return 0;
    }
  } catch (po::error &e) {
    std::cerr << e.what() << std::endl;
    std::cout << "USAGE: " << av[0] << std::endl;
    std::cout << description << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<std::string> in_files = vm["in_files"].as<std::vector<std::string> >();
  std::string out_file = vm["out_file"].as<std::string>();

  std::ofstream ofile;
  std::ifstream ifile;
  int N = 0;
  int i = 0;
  double eps, rho;
  std::vector<double> rho_vals;
  std::vector<double> eps_vals;

  ifile.open(in_files[0].c_str());
  if (ifile.is_open()) {
    N = std::count(std::istreambuf_iterator<char>(ifile),
               std::istreambuf_iterator<char>(), '\n');
    ifile.close();
  }

  std::cout << "number of points: " << N << std::endl;

  rho_vals.resize(N);
  eps_vals.resize(N);

  for (auto file : in_files) {

    std::cout << "opening file: " << file << std::endl;

    ifile.open(file);
    if (ifile.is_open()) {
      while (ifile >> eps >> rho) {
        rho_vals[i] += rho;
        eps_vals[i] = eps;
        i++;
      }
      ifile.close();
      i = 0;
    }

  }


  ofile.open(out_file);
  if (ofile.is_open()){
    for (int i = 0; i < rho_vals.size(); i++){
      rho = rho_vals[i];
      rho /= in_files.size();

      ofile << eps_vals[i] << " " << rho << std::endl;
    }
    ofile.close();
  }


  return 0;
}
