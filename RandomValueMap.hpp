#ifndef RANDOMVALUEMAP
#define RANDOMVALUEMAP


class RandomValueMap {
public:
	double operator()(int x);
	void clear();
private:
	std::unordered_map<int,double> rand_vec;	


};

#endif
