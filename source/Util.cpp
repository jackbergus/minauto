#include"../include/minauto/Util.hpp"

bool global_setting::html = false;
int global_setting::verbosity = 1;
bool global_setting::interval = false;
bool global_setting::check_minimisation = false;
bool global_setting::prune = false;


long double relative_precision = 1E-6;
long double absolute_precision = std::numeric_limits<long double>::epsilon();

double getRandomDouble(unsigned long long K) {
  double result( (rand() % (K + 1) ) / ((double) K));
  return result;
}

std::string toString(float f) {
	static char buffer[33];
	snprintf(buffer,33,"%G",f);
	return buffer;	
}

std::string toString(double f) {
	static char buffer[33];
	snprintf(buffer,33,"%lG",f);
	return buffer;	
}

std::string toString(long double f) {
	static char buffer[33];
	snprintf(buffer,33,"%LG",f);
	return buffer;	
}

std::string toString(int i) {
  static char buffer[33];
  snprintf(buffer,33,"%d",i);
  return std::string(buffer);
}

std::string toString(unsigned int i) {
  static char buffer[33];
  snprintf(buffer,33,"%d",i);
  return std::string(buffer);
}

std::string toString(mpfr::mpreal f) {
  return f.toString();
}

std::string toString(const mpfr::mpreal& f) {
  return f.toString();
}

void progressBar(int j, int n) {
	float ratio = j / (float) n;
	const int bar_length = 20;
	
	if(global_setting::verbosity >= 1) {
		std::cout << "  |";
	
		for(int i=0; i<bar_length; ++i) {
			if(i < ratio * bar_length)
				std::cout << "=";
			else
				std::cout << " ";
		}
		std::cout << "|\r" << std::flush;
	}
}

