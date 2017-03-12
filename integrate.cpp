#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <boost/program_options.hpp>

int main(int ac, char** av) {

  namespace po = boost::program_options;
  po::options_description description("Options");
  description.add_options()
    ("in_file,i", po::value<std::string>()->required(), "input file");
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

  std::string in_file = vm["in_file"].as<std::string>();
  std::vector<std::pair<double, double> > rho0;

  std::ifstream dosfile;
  dosfile.open(in_file);
  {
    double w,rho;
    while(dosfile >> w >> rho)
      rho0.push_back(std::make_pair(w,rho));
  }
  dosfile.close();

  std::cout << "rho0 loaded from file: " << rho0.size() << "points." << std::endl;

  double preve = rho0.front().first;
  double prevf = rho0.front().second;
  double integral = 0;
  for(auto const& rho0pt : rho0)
    {
      double e = rho0pt.first;
      double f = rho0pt.second;
      integral += 0.5*(prevf + f)*(e - preve);
      preve = e;
      prevf = f;
    }
  std::cout << std::setprecision(16) << "Density of states integrates to:" << integral << std::endl;

  return 0;
}
