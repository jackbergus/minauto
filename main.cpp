#include "Util.hpp"
#include "Timer.hpp"

#include "BLAS.hpp"
#include "Automaton.hpp"

//#include "GramSchmidt.hpp"
#include "WordTree.hpp"
#include "Worklist.hpp"
#include "RandomValueMap.hpp"

#include "GramSchmidt.hpp"
#include "Minimisation.hpp"
#include "EquivAlgorithms.hpp"
#include "Inflate.hpp"
#include "Mutate.hpp"

#include <ctime>




enum Mission {
	  minimisation,
	  equivalence
};




std::string command_line_options =
  "Minimiser and equivalence checker for probabilistic automata\n"\
  "University of Oxford 2011-2013\n"\
  "\n"\
  "Manual\n"\
  "\n"\
  "equiv <options>* <automaton1> <automaton2>\n"\
  "\n"\
  "options: \n"\
  "\n"\
  "------ Export --------\n"\
  "  - matlab                    ... export automaton to Matlab\n"\
  "------ Precision -----\n"\
  "  -onion                      ... continue with Gram-Schmidt-projected vector\n"\
  "  -absolute_precision <value> ... absolute numerical precision\n"\
  "  -relative_precision <value> ... relative numerical precision\n"\
  "  -normalise                  ... normalise in Gram Schmidt\n"\
  "------ Minimisation & Reduction ----\n"\
  "  -backward_reduction         ... backward reduction only\n"\
  "  -forward_reduction          ... forward reduction only\n"\
  "  -minimise               ... backward & forward reduction"\
  "------ Equivalence ----\n"\
  "  -ce                         ... switch on counterexample generation and checking\n"\
  "  -forward                    ... Tzeng's forward algorithm\n"\
  "  -backward                   ... Henzinger and Raskin's backward algorithm\n"\
  "  -randomised_forward         ... randomised forward algorithm \n"\
  "  -randomised_backward        ... randomised backward algorithm \n"\
  "  -html                       ... produce HTML output only \n"\
  "------ Mutation ----\n"\
  "  -forward_inflate  <steps>   ... forward-inflate the automaton with additional state\n"\
  "  -backward_inflate <steps>   ... backward-inflate the automaton with additional state\n"\
  "  -mutate <steps>             ... mutate automaton to inequivalent similar automaton\n"\
;
  

struct ModStep {

  enum Operation { forward, backward, mutate} operation;
  unsigned steps;

  ModStep(Operation o, unsigned s) : operation(o) , steps(s) {}
};

double check_absolute_precision = 1E-06,
       check_relative_precision = 1E-06;

enum OutputFormat { Matlab, aiSee, Generic } output_format;

std::list<ModStep> mod_sequence;


template<typename Automaton>
bool check(	MinimisationSettings& msettings,Settings& settings, Automaton& a, Automaton& b, Precision<typename Automaton::Weight>& prec) {

  /*
	typedef ::Automaton<mpfr::mpreal> MPAutomaton ;

	  mpfr::mpreal::set_default_prec(256);

	  mpfr::mpreal absolute_precision =  (check_absolute_precision == 1e-6) ? std::numeric_limits<mpfr::mpreal>::epsilon() : check_absolute_precision; 

	  Precision<mpfr::mpreal> prec(absolute_precision,check_relative_precision);
  */
	Automaton c(a);
	Automaton d(b);

	Automaton diff(c,d);
	  
	unsigned k = a.getNrOfStates();

	mpfr::mpreal maximum = kmax<Automaton,true> ( diff, k, settings, prec);

	if(global_setting::verbosity >= 1) {
	std::cout << "Maximum difference " << maximum << " after " << k << " steps " <<std::endl;
	}
	return true; // equiv(c,d,settings,prec);
}


template<typename Automaton>
int run(Mission& mission,
	std::vector<std::string>& automata_files,
	MinimisationSettings& msettings,
	Settings& settings,
	int precision, // how many bits
	int inflation_precision, // how many bits during inflation
	Precision<typename Automaton::Weight>& prec,
	std::list<ModStep>& mod_sequence,
	OutputFormat output_format,
	bool print) {
	unsigned auto_counter = automata_files.size();

	Automaton a;
	Automaton b;

	switch(mission) {
	  case minimisation: {

      std::cout << "Minimiser" << std::endl;

		  if(auto_counter>2) {
		  		  std::cout << "Minimisation expects one or two automata as an argument" << std::endl;
		  		  std::cout << " " << auto_counter << " arguments were provided" << std::endl;
		  }

		  a.parse(automata_files[0].c_str());
	      a.computeMatrices();

		  if(auto_counter == 2) {
			  b.parse(automata_files[1].c_str());
			  b.computeMatrices();
		  }

	      if(inflation_precision!=-1) {
	      		  mpfr::mpreal::set_default_prec(inflation_precision);
	      }

	      foreach(ModStep m, mod_sequence) {
	      	    switch(m.operation) {
	      			case ModStep::forward:
	      				forwardInflate(a,m.steps);
	      				break;
	      			case ModStep::backward:
	      				backwardInflate(a,m.steps);
	      				break;
	      			case ModStep::mutate:
	      				mutateAutomaton(a,m.steps);
	      				break;
	      	    }
	      }

	      if(inflation_precision!=-1) {
	      		  mpfr::mpreal::set_default_prec(precision);
	       }

	      if(print)
	      switch(output_format)  {
			case aiSee: {

				std::string auto1_file("auto.gdl");

				std::ofstream file1 (auto1_file.c_str());

				a.aiSee(file1);
				exit(0);
				}
			case Matlab: {
				a.writeMatlab("auto.m");
				exit(0);
			}
			case Generic: {
			  a.write("auto.aut");
				  exit(0);
			}
	      }

	      switch(auto_counter) {
			  case 1: 
			    minimise(a,b,msettings,prec);
			    std::cout<<"   # states: " << b.getNrOfStates() 
			             << " # trans: "   << b.getNrOfNonzeros() 
			             << " # letter " << alphabet.size()<< std::endl;
					  b.write("mini.aut");

					  if(global_setting::check_minimisation) {
					    check<Automaton>(msettings,settings,a,b,prec);
				          }

			  break;
			  case 2: {
				  Automaton c(a,b);
				  c.write("minus.aut");

				  Automaton d;
				  minimise(c,d,msettings,prec);
				  std::cout<<"   # states: " << d.getNrOfStates() << " # trans: " << d.getNrOfNonzeros() << " # letter " << alphabet.size()<< std::endl;
				  d.write("mini.aut");

				  if(global_setting::check_minimisation) {
				    check<Automaton>(msettings,settings,c,d,prec);
				  }


			  }
	      }

	      if(global_setting::verbosity >= 1) {
	    	  std::cout << "absolute precision " << prec.absolute_precision << std::endl;
	      }
	     break;
	  }
	  case equivalence:  {


		  if(auto_counter!=2) {
			  std::cout << "Equivalence check expects exactly two automata as an argument" << std::endl;
			  std::cout << " " << auto_counter << " arguments were provided" << std::endl;
		  }

	      a.parse(automata_files[0].c_str());
	      a.computeMatrices();
	      b.parse(automata_files[1].c_str());
	      b.computeMatrices();

	      if(inflation_precision!=-1) {
			  mpfr::mpreal::set_default_prec(inflation_precision);
	      }

	      foreach(ModStep m, mod_sequence) {
		    switch(m.operation) {
				case ModStep::forward:
					forwardInflate(a,m.steps);
					break;
				case ModStep::backward:
					backwardInflate(a,m.steps);
					break;
				case ModStep::mutate:
					mutateAutomaton(a,m.steps);
					break;
		    }
	      }

	      if(inflation_precision!=-1) {
	      		  mpfr::mpreal::set_default_prec(precision);
	       }

	      if(global_setting::verbosity >= 1) {
	      	//printf("Memory (automata): %ld bytes\n",a.getMem() + b.getMem());
	      }
	      return equiv(a,b,settings,prec);
	    }
	    default:
	      std::cout << "Not clear what to do" << std::endl;
	      break;
	    
  }
  return 0;
}


int main(int argc, char* argv[]) {

  double absolute_precision = 1E-16,
		 relative_precision = 1E-06;



	
  enum RealType {
	  float32,
	  float64,
	  float80,
	  mp,
	  rational
  } reals = float32;
	

  Mission mission;



  std::vector<std::string> automata_files;
  
  MinimisationSettings msettings;
  Settings settings;

  msettings.lz = false;
  msettings.direction = MinimisationSettings::both;
  msettings.pivot = false;
  msettings.reortho = false;


  settings.pivot = false;
  settings.reortho = false;


  settings.ce = false;
  settings.onion = false;
  settings.normalise = false;
  global_setting::verbosity = 1;
  global_setting::prune = false;


  bool print = false;
  bool custom_absolute_precision = false;

  int precision = 128;
  int inflation_precision = -1;
  
  /* initialize random seed: */
  srand ( 41 );

  if(argc == 1) {
    std::cerr <<command_line_options<<std::endl;
    exit(1);
  }

  for(int i=1; i<argc; ++i) {
    const std::string arg = argv[i];
    if(arg == "-ce") {
      settings.ce = true;
    }
    else if(arg == "-classical") {
          msettings.classical=true;
    }
    else if(arg == "-minimise") {
          mission = minimisation;
    }
    else if(arg == "-equivalence") {
              mission = equivalence;
    }
    else if(arg == "-pivot") {
              msettings.pivot = settings.pivot = true;
    }
    else if(arg == "-reortho") {
              msettings.reortho = settings.reortho = true;
    }
    else if (arg == "-check") {
    	global_setting::check_minimisation = true;
    }
    else if (arg == "-prune") {
        	global_setting::prune = true;
    } else if ( arg == "-html" ) {
      global_setting::html = true;
      global_setting::verbosity = 0;
    } else if(arg == "-float32") {
      reals = float32;
    } else if(arg == "-float64") {
      reals = float64;
    } else if(arg == "-float80") {
        reals = float80;
  } else if(arg == "-mp") {
          reals = mp;
        }

  else if(arg == "-interval") {
	global_setting::interval = true;

  } else if(arg == "-forward") {
      settings.method = Settings::forward;
      msettings.direction = MinimisationSettings::forward;
	}
    else if(arg == "-lz") {
              msettings.lz = true;
        }
    else if(arg == "-householder") {
          msettings.method = MinimisationSettings::householder;
    }
    else if( arg == "-backward") {
      settings.method = Settings::backward;
      msettings.direction = MinimisationSettings::backward;
    } 
    else if( arg == "-randomised_backward") {
      settings.method = Settings::randomised_backward;
   }
    else if( arg == "-randomised_forward") {
      settings.method = Settings::randomised_forward;
    } 
    else if ( arg == "-precision" ) {
      if( i < argc - 1 && sscanf(argv[i+1],"%d",&precision) == 1) {
    	  ++i;
      }
    }
    else if ( arg == "-inflation_precision" ) {
          if( i < argc - 1 && sscanf(argv[i+1],"%d",&inflation_precision) == 1) {
        	  ++i;
          }
    }
    else if ( arg == "-absolute_precision" ) {
      custom_absolute_precision = true;
      if( i < argc - 1 && sscanf(argv[i+1],"%lG",&absolute_precision) == 1) {
    	  ++i;
      }
    } else if ( arg == "-relative_precision" ) {
      if( i < argc - 1 && sscanf(argv[i+1],"%lG",&relative_precision) == 1) {
    	  ++i;
      }
    }
    else if ( arg == "-check_absolute_precision" ) {
      if( i < argc - 1 && sscanf(argv[i+1],"%lG",&check_absolute_precision) == 1) {
    	  ++i;
      }
    } else if ( arg == "-check_relative_precision" ) {
      if( i < argc - 1 && sscanf(argv[i+1],"%lG",&check_relative_precision) == 1) {
    	  ++i;
      }
    }

    else if ( arg == "-verbosity" ) {
      if( i < argc - 1 && sscanf(argv[i+1],"%d",&global_setting::verbosity) == 1) {
    	  ++i;
      }
    }
    else if ( arg == "-forward_inflate" ) {
      int steps;
      if( i < argc - 1 && sscanf(argv[i+1],"%d",&steps) == 1) {
    	  ++i;
    	  ModStep i_step(ModStep::forward,steps);
    	  mod_sequence.push_back(i_step);
      }
    }
    else if ( arg == "-backward_inflate" ) {
      int steps;
      if( i < argc - 1 && sscanf(argv[i+1],"%d",&steps) == 1) {
    	  ++i;
    	  ModStep i_step(ModStep::backward,steps);
    	  mod_sequence.push_back(i_step);
      }
    }
    else if ( arg == "-mutate" ) {
      int steps;
      if( i < argc - 1 && sscanf(argv[i+1],"%d",&steps) == 1) {
    	  ++i;
    	  ModStep i_step(ModStep::mutate,steps);
    	  mod_sequence.push_back(i_step);
      }
    } 
    else if( arg == "-onion") {
      settings.onion = true; 
    }
    else if( arg == "-normalise") {
      settings.normalise = true;
    } 
    else if ( arg == "-aisee") {
		print = true;
		output_format = aiSee;
    }	
    else if ( arg == "-print") {
		print = true;
		output_format = Generic;
    }	
    else if ( arg == "-matlab") {
		print = true;
		output_format = Matlab;
    }	
    else {
      automata_files.push_back(arg.c_str());
    }
  }

  
  switch(reals) {
  case float32:
	  {
		  absolute_precision =  (custom_absolute_precision == false) ? std::numeric_limits<float>::epsilon() : absolute_precision; 
		  Precision<float> prec(absolute_precision,relative_precision);
		  return run<Automaton<float> >(mission,automata_files,msettings,settings,precision,inflation_precision,prec,mod_sequence,output_format, print);
	  }
	  break;
	  
  case float64:
	  {
		  absolute_precision =  (custom_absolute_precision == false) ? std::numeric_limits<double>::epsilon() : absolute_precision; 
		  Precision<double> prec(absolute_precision,relative_precision);
		  return run<Automaton<double> >(mission,automata_files,msettings,settings,precision,inflation_precision,prec,mod_sequence,output_format, print);
	  }
	  break;
  case float80:
	  {
		  absolute_precision =  (custom_absolute_precision == false) ? std::numeric_limits<long double>::epsilon() : absolute_precision; 
		  Precision<long double> prec(absolute_precision,relative_precision);
		  return run<Automaton<long double> >(mission,automata_files,msettings,settings,precision,inflation_precision,prec,mod_sequence,output_format, print);
	  }
	  break;
  case mp:
	  {
		  mpfr::mpreal::set_default_prec(precision);
		  mpfr::mpreal ab(absolute_precision);
		  Precision<mpfr::mpreal> prec(((custom_absolute_precision == false) ? mpfr::machine_epsilon() : ab),relative_precision);
		  return run<Automaton<mpfr::mpreal> >(mission,automata_files,msettings,settings,precision,inflation_precision,prec,mod_sequence,output_format, print);
	  }
	  break;
  default:
	  break;
  }
  return 0;
  
}
