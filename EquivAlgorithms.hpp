#ifndef EQUIVALGORITHMS
#define EQUIVALGORITHMS

struct Settings {

  Settings () : method (forward), full (false) {}

  enum Method {
    forward,
    backward,
    randomised_backward,
    randomised_forward
  } method;

  bool ce;
  bool onion;  // continue with projected or original vector - onion=true may improve numerical stability 
  bool normalise;  

  bool full;

  bool pivot;
  bool reortho;
};

template <typename Automaton>
struct Result {
  typedef typename Automaton::Weight Weight;
	
	
  bool equivalent;  // automata are equivalent or not
  Word ce;          // counterexample if automata are inequivalent (and CE generation is switched on)
  Weight a;         // probability of counterexample word in automaton A
  Weight b;         // probability of counterexample word in automaton B
  int node;         // node in the language tree with maximal counterexample


  /* statistics */
  int dimensions;
  long iterations;
  long time;                           // time in seconds
};


/** Equivalence checking based on Tzeng's algorithm
    template parameters
    - Automaton ... type of automaton which determines the underlying data structures and numeric types
    - forward   ... direction of exploration
    function parameters
    - a        ...  automaton
    - b        ...  automaton
    - settings ...  see above
    - prec     ...  arithmetic precision and threshold values
    - r        ...  result of the equivalence check
 */
template <typename Automaton, bool forward>
bool Tzeng ( Automaton& a, Automaton& b,
             const Settings& settings, 
             Precision<typename Automaton::Weight>& prec, Result<Automaton>& r) {

  if(global_setting::verbosity >= 1) {
  	std::cout << (forward ? "Forward" : "Backward" ) << " Tzeng " <<std::endl;
  }

  typedef typename Automaton::Weight Scalar;
  typedef typename Automaton::Vector Vector;
  typedef typename Automaton::Matrix Matrix;
  typedef GramSchmidt<Automaton> GS;
  
  bool result = true;

  long inner_loop_counter = 0;

  using namespace boost::numeric::ublas;
  typedef typename Vector::const_iterator VectorConstIterator;
  typedef boost::numeric::ublas::vector_range<Vector> VectorRange;

  if( norm_1(a.getStart(forward)) == 0.0 && norm_1(a.getStart(forward)) == 0.0) {
    // Both automata accept the empty language
    return true;
  }

  // check for the empty word
  Scalar proda (inner_prod(a.initial,a.accepting));
  Scalar prodb (inner_prod(b.initial,b.accepting));
  
  if( !prec.isEqual(proda,prodb) ) {
    // Empty word has different probability
    return false;
  }

  // start worklist algorithm
  
  Worklist<DWPair<Vector> > worklist;
  unsigned na = a.getStart(forward).size();
  unsigned nb = b.getStart(forward).size();
  unsigned n = na + nb;

  GS gs(na + nb, prec, settings.normalise,settings.pivot,settings.reortho);

  Vector start (na + nb);
  subrange(start,0,na) = a.getStart(forward);
  subrange(start,na,na+nb) = b.getStart(forward);
  start *= 1/norm_2(start); 

  gs.addToBasis(start);


  WordTree tree;

  DWPair<Vector> dwp (a.getStart(forward), b.getStart(forward), tree.addNode(0,-1)) ;

  worklist.add(dwp);

  Vector q (na + nb);
  q.clear();

  unsigned int l = 0; // records the unique index of the current work list item

  Scalar max_diff(0);
  
  while (worklist.size() > 0) {
    dwp = worklist.top();
    const Vector va = dwp.ua;
    const Vector vb = dwp.ub;
    worklist.pop();

    // get the enabled letters
    std::unordered_set<int> letters;

    a.getLetters(va,letters,forward);
    b.getLetters(vb,letters,forward);

    progressBar(l,n);

    foreach(int sigma, letters) {
      Vector ua ( na ) ;
      ua.clear();
      a(sigma,ua,va,forward);

      Scalar pa(prec_inner_prod(ua,a.getFinish(forward)));

      Vector ub ( nb );
      ub.clear();
      b(sigma,ub,vb,forward);

      Scalar pb(prec_inner_prod(ub,b.getFinish(forward)));
     
      if( prec.isEqual(pa,pb)) {

		subrange(q,0,na) = ua;
		subrange(q,na,na+nb) = ub;


		++inner_loop_counter;

		gs.project(q);


		if(gs.addToBasis(q) && l < n) {
		  DWPair<Vector> dwp_new;
		  dwp_new.ua = settings.onion ? subrange(q,0,na) : ua;
		  dwp_new.ub = settings.onion ? subrange(q,na,na+nb) : ub;

		  if(settings.ce) { // counterexample generation enabled?
			dwp_new.node = tree.addNode(sigma,dwp.node);
		  }

		  ++l;

		  worklist.add(dwp_new);
		}
	} 	else {
		result = false;

	    if(global_setting::verbosity >= 1) {
			std::cout << "Difference " << fabs(pb - pa) << std::endl;
		}
		
		int node = 0;
		
		if(settings.ce) { // counterexample generation enabled?
			node = tree.addNode(sigma,dwp.node);
			
			Scalar diff(fabs(pa - pb));
			
			if(diff > max_diff) {
			    max_diff = diff;
			    r.node = node;
			}
			
		}
		
		
		if(!settings.full) {
		  worklist.clear(); 
		}
		else {
		  subrange(q,0,na) = ua;
		  subrange(q,na,na+nb) = ub;


		  ++inner_loop_counter;

		  gs.project(q);


		  if(gs.addToBasis(q) && l < n) {
		    DWPair<Vector> dwp_new;
		    dwp_new.ua = settings.onion ? subrange(q,0,na) : ua;
		    dwp_new.ub = settings.onion ? subrange(q,na,na+nb) : ub;

		    if(settings.ce) { // counterexample generation enabled?
			  dwp_new.node = node;
		    }

		    ++l;

		    worklist.add(dwp_new);
		  }
		}
		break;
      }
    }
    ++l;
  }

  if(result == false && settings.ce) { // counterexample generation enabled?
  	tree.getWord(r.node,r.ce,forward);		    	  
  }

  // store the result
  r.equivalent  = result;
  r.dimensions  = gs.getDimension();
  r.iterations  = inner_loop_counter;

  return result;
}

template <typename Automaton, bool forward>
typename Automaton::Weight kmax ( 
	Automaton& a, 
	unsigned k, 
	const Settings& settings, 
        Precision<typename Automaton::Weight>& prec) {

  if(global_setting::verbosity >= 1) {
  	std::cout << (forward ? "Forward" : "Backward" ) << " Exploration " <<std::endl;
  }

  typedef typename Automaton::Weight Scalar;
  typedef typename Automaton::Vector Vector;
  typedef typename Automaton::Matrix Matrix;
  
  using namespace boost::numeric::ublas;
  typedef typename Vector::const_iterator VectorConstIterator;
  typedef boost::numeric::ublas::vector_range<Vector> VectorRange;

  if( norm_1(a.getStart(forward)) == 0.0 ) {
    // zero automaton
    return 0;
  }

  // check for the empty word?
  
  // start worklist algorithm
  
  Worklist<DW<Vector> > worklist;
  unsigned n = a.getStart(forward).size();

  Vector start ( a.getStart(forward) );

  WordTree tree;

  DW<Vector> dw(start, tree.addNode(0,-1), 0);


  worklist.add(dw);

  Vector q (n);
  q.clear();

  unsigned int l = 0; // records the unique index of the current work list item

  Scalar maximum (0);
  
  unsigned current_depth = 0;
  int maximal_node = -1;

  unsigned nr_mult=0;
  unsigned max_depth=0;


  while (worklist.size() > 0) {

    dw = worklist.top();
    current_depth = dw.depth;

    worklist.pop();

    if(current_depth >= k)
      continue;

    const Vector v = dw.u;


    // get the enabled letters
    std::unordered_set<int> letters;

    a.getLetters(v,letters,forward);


    progressBar(current_depth,k);

    foreach(int sigma, letters) {
      Vector u ( n ) ;
      a(sigma,u,v,forward);

      ++nr_mult;

      Scalar p(fabs(prec_inner_prod(u,a.getFinish(forward))));

      if(p > maximum && p != std::numeric_limits<double>::infinity()) {
      	maximum = p;
      	max_depth=dw.depth;
		if(global_setting::verbosity >= 1) {
      		std::cout << "Curr max " << p << std::endl;
      	}
      }

      int node = 0;

      if(settings.ce) { // counterexample generation enabled?
	      node = tree.addNode(sigma,dw.node);
      }

      if(!prec.isZero(norm_inf(u))) {
        dw.node = node;
        dw.u = u;
        dw.depth = current_depth + 1;
     
	      worklist.add(dw);
      }     
    }
    ++l;
  }
  
  if(global_setting::verbosity >= 1) {
  	std::cout << "kmax nr of steps " << nr_mult << " max at depth " << max_depth << std::endl;
  }

  Word word;

  if(settings.ce) { // counterexample generation enabled?
  	tree.getWord(maximal_node, word,forward);		    	  
  }

  return maximum;
}





/*
 * compute image under automaton with respect to given random vector 
 * see also: Randomised2
 */
template <typename Automaton>
bool computeRandomisedImage(const typename Automaton::Vector& x, Automaton& a, RandomValueMap& R, typename Automaton::Vector& output) {
	typedef typename Automaton::Vector Vector;
	// get the enabled letters
	unsigned na = x.size();
    std::unordered_set<int> enabled;
    a.get_enabled(x,enabled);

    for(std::unordered_set<int>::const_iterator it=enabled.begin(); it!=enabled.end();++it) {
	  int sigma(*it);
	  Vector vsigma(na);
	  a(sigma,vsigma,x,true);
	  vsigma *= R(sigma);
	  output +=vsigma;
    }
    return true;
}

/*
template <class T>
boost::shared_future<T> submitJob(boost::function<T (void)> const& fn) {
	boost::promise<T> prom;
	prom.set(fn());
	return boost::shared_future<T>(prom);
}	



struct Worker {

	boost::unique_future<void> work(boost::function<void ()> fn) {
		boost::packaged_task<void> pt(fn);
		task = boost::move(pt);
		condition.notify_all();
		return task.get_future();
	}

	void run() {
		boost::mutex::scoped_lock lock(mutex);
		while(true) {

			std::cout<<"Worker running and waiting for condition" << std::endl;
			condition.wait(lock);
			std::cout<<"Condition fulfilled" << std::endl;

			task();
		}
	}
	boost::packaged_task<void> task;
	boost::mutex mutex;
	boost::condition_variable condition;
};

*/ 

template <typename Automaton>
bool Randomised2 ( Automaton& a, Automaton& b, const Settings& settings, Precision<typename Automaton::Weight>& prec, Result<Automaton>& r) {
  typedef typename Automaton::Weight Scalar;
  typedef typename Automaton::Vector Vector;
	
  bool result = true;

  long inner_loop_counter = 0;
  unsigned terminating_iteration = 0;
  


  using namespace boost::numeric::ublas;
  typedef typename Vector::const_iterator VectorConstIterator;
  typedef boost::numeric::ublas::vector_range<Vector> VectorRange;


  if( norm_1(a.accepting) == 0.0 && norm_1(b.accepting) == 0.0) {
    return true;
  }

  // check for the empty word
  Scalar proda (prec_inner_prod(a.initial,a.accepting));
  Scalar prodb (prec_inner_prod(b.initial,b.accepting));
  
  if( !prec.isEqual(proda,prodb) ) {
    return false;
  }

  // start worklist algorithm
  unsigned na = a.initial.size();
  unsigned nb = b.initial.size();
  unsigned n = na + nb;

  Vector va (a.initial), vb(b.initial);
  RandomValueMap R;


  std::vector<Vector> ma, mb;

  ma.push_back(va);
  mb.push_back(vb);

  if(global_setting::verbosity >= 1) {
  	std::cout << "Maximal nr of iterations: " << n << std::endl;
  }  
  //boost::thread_group threads;

  /*
  Worker worker1, worker2;

  threads.create_thread(boost::bind( &Worker::run, &worker1 ));
  threads.create_thread(boost::bind( &Worker::run, &worker2 ));
  

  worker1.work(boost::bind( &f,"Making worker 1 busy\n"));
  worker2.work(boost::bind( &f,"Making worker 2 busy\n"));
  worker1.work(boost::bind( &f,"Making worker 1 busy\n"));
  worker2.work(boost::bind( &f,"Making worker 2 busy\n"));
   */

  for(unsigned i=0; result && i<n; ++i) {
    Vector va_new (na);
    va_new.clear();
    Vector vb_new (nb);
    vb_new.clear();

    computeRandomisedImage<Automaton>(va, a, R, va_new);
    va = va_new;
    Scalar pa(prec_inner_prod(va,a.accepting));

    computeRandomisedImage<Automaton>(vb, b, R, vb_new);
    vb = vb_new;

    Scalar pb(prec_inner_prod(vb,b.accepting));

    R.clear();


    if(settings.ce) {
      // remember vectors
      ma.push_back(va);
      mb.push_back(vb);
    }

    inner_loop_counter+=2;
      

    if(!prec.isEqual(pa,pb)) {
	  if(global_setting::verbosity >= 1) {
      	std::cout << "Difference " << fabs(pa - pb) << std::endl;
      }

      terminating_iteration = i;
      result = false;
    }
  }

  if(global_setting::verbosity >= 1) {
  	std::cout << "Waiting for threads to terminate" << std::endl;
  }
  
  //threads.join_all();

  if(global_setting::verbosity >= 1) {
  if(result == false)
  printf("%d iterations\n",terminating_iteration+1);
  }
  
  // analyse the counterexample
  if(result == false && settings.ce) {

    Vector ua(a.accepting), ub(b.accepting);

    if(terminating_iteration + 2 != ma.size()) {
    	if(global_setting::verbosity >= 1) {
		  std::cout << "Failed to generate CE" << std::endl;
		}
    }

	  if(global_setting::verbosity >= 1) {
	    std::cout << "A difference was found at iteration " << terminating_iteration 
	              <<". Try to generate a counterexample word." << std::endl;
	  }
	  
    for(int i= terminating_iteration; i>-1; --i) {
      // find an action where there's a difference

      // get the enabled letters
      std::unordered_set<int> incoming;
      a.get_incoming(ua,incoming);
      b.get_incoming(ub,incoming);

      Scalar discrepancy = 0;
      int letter = -1;

      Vector uat (na);
      Vector ubt (nb);

      // find the action with the maximal relative discrepancy
      for(std::unordered_set<int>::const_iterator it=incoming.begin(); it!=incoming.end();++it) {
		int sigma(*it);

	        uat.clear();
		a(sigma,uat,ua,false);
	
  		ubt.clear();
		b(sigma,ubt,ub,false);

		Scalar pa = prec_inner_prod(uat,ma[i]);
		Scalar pb = prec_inner_prod(ubt,mb[i]);

	        Scalar diff(fabs(pa - pb));

		if(discrepancy < diff) {
			discrepancy = diff;
		  	letter = sigma;
		}	
      }



      if(letter == -1) {
	      if(global_setting::verbosity >= 1) {
			std::cout << "Randomised: Error in CEX analysis steps "<< i << std::endl;
		  }
		break;
      } else {
	      r.ce.push_front(letter);
              ua = uat;
              ub = ubt;
      }

    }
  }

  // store the result
  r.equivalent  = result; 
  r.dimensions  = terminating_iteration+1;
  r.iterations  = terminating_iteration+1;

  return result;
}

template <typename Automaton>
bool Randomised ( Automaton& a, Automaton& b, const Settings& settings, Precision<typename Automaton::Weight>& prec, Result<Automaton>& r) {
  typedef typename Automaton::Weight Scalar;
  typedef typename Automaton::Vector Vector;
	  
  bool result = true;
  long inner_loop_counter = 0;
  unsigned terminating_iteration = 0;

  using namespace boost::numeric::ublas;
  typedef typename Vector::const_iterator VectorConstIterator;
  typedef boost::numeric::ublas::vector_range<Vector> VectorRange;


  if( norm_1(a.accepting) == 0.0 && norm_1(b.accepting) == 0.0) {
	  if(global_setting::verbosity >= 1) {
	    std::cout<< "Both automata accept the empty language"<<std::endl;
	  }
    return true;
  }

  // check for the empty word
  Scalar proda (prec_inner_prod(a.initial,a.accepting));
  Scalar prodb (prec_inner_prod(b.initial,b.accepting));
  


  if( !prec.isEqual(proda,prodb) ) {
    if(global_setting::verbosity >= 1) {
	    std::cout<< "Empty word has different probability"<<std::endl;
	}
    return false;
  }

  // start worklist algorithm
  unsigned na = a.accepting.size();
  unsigned nb = b.accepting.size();
  unsigned n = na + nb;

  Vector va (a.accepting), vb(b.accepting);
  RandomValueMap R;

  if(global_setting::verbosity >= 1) {
	  std::cout << "Maximal nr of iterations: " << n << std::endl;
  }

  std::vector<Vector> ma, mb;

  ma.push_back(va);
  mb.push_back(vb);

  for(unsigned i=0; result && i<n; ++i) {
    Vector va_new (na);
    Vector vb_new (nb);

    // get the incoming letters
    std::unordered_set<int> incoming;

    a.get_incoming(va,incoming);

    for(std::unordered_set<int>::const_iterator it=incoming.begin(); it!=incoming.end();++it) {
	  int sigma(*it);
	  Vector vsigma (na);
	  a(sigma,vsigma,va,false);
	  vsigma *= R(sigma);
	  va_new +=vsigma;
    }
    va = va_new;

    Scalar pa(prec_inner_prod(va,a.initial));

    incoming.clear();
    b.get_incoming(vb,incoming);
 

    for(std::unordered_set<int>::const_iterator it=incoming.begin(); it!=incoming.end();++it) {
	  int sigma(*it);
	  Vector vsigma (nb);
	  b(sigma,vsigma,vb,false);
	  vsigma *= R(sigma);
	  vb_new +=vsigma;
    } 
    vb = vb_new;

    Scalar pb(prec_inner_prod(vb,b.initial));

    R.clear();

    if(settings.ce) {
      // remember vectors
      ma.push_back(va);
      mb.push_back(vb);
    }

    inner_loop_counter+=2;

    if(!prec.isEqual(pa,pb)) {
	  if(global_setting::verbosity >= 1) {
	      std::cout << "Difference " << fabs(pa - pb) << std::endl;
	  }

      terminating_iteration = i;
      result = false;
    }
  }

  // get the counterexample
  if(result == false && settings.ce) {
	  if(global_setting::verbosity >= 1) {
	    std::cout << "A difference was found at iteration " 
	    		  << terminating_iteration <<". Try to generate a counterexample word." 
	    		  << std::endl;
	  }
    Vector ua(a.initial), ub(b.initial);

    if(terminating_iteration + 2 != ma.size()) {
  		if(global_setting::verbosity >= 1) {
	  		std::cout << "ma: Whisky tango foxtrot?" << std::endl;
	  	}
    }

    for(int i= terminating_iteration; i>-1; --i) {



      // find an action where there's a difference

      // get the enabled letters
      std::unordered_set<int> enabled;
      a.get_enabled(ua,enabled);
      b.get_enabled(ub,enabled);

      Scalar discrepancy = 0;
      int letter = -1;
      Vector uat(na);
      Vector ubt(nb);

      for(std::unordered_set<int>::const_iterator it=enabled.begin(); it!=enabled.end();++it) {
    	  int sigma(*it);
	

	  uat.clear();
	  a(sigma,uat,ua,true);

	  ubt.clear();
	  b(sigma,ubt,ub,true);
		
	  Scalar pa = prec_inner_prod(uat,ma[i]);
	  Scalar pb = prec_inner_prod(ubt,mb[i]);

          Scalar diff(fabs(pa - pb));

	  if(discrepancy < diff) {
	  	discrepancy = diff;
	  	letter = sigma;
	   }
      }

      if(letter == -1) {
		  if(global_setting::verbosity >= 1) {
			std::cout << "Randomised: Error in CEX analysis steps "<< i << std::endl;
		  }
		break;
      } else {
	      r.ce.push_back(letter);
              ua = uat;
              ub = ubt;
      }
    }
  }

  if(global_setting::verbosity >= 1) {
	  if(result == false)
	  printf("%d iterations\n",terminating_iteration+1);
  }

  // store the result
  r.equivalent  = result; 
  r.dimensions  = terminating_iteration+1;
  r.iterations  = terminating_iteration+1;

  return result;
}

template<typename Automaton>
void createHTMLReport(const Result<Automaton>& result, std::string& report) {

	
	report += "<!-- Probabilistic Equivalence Checker (University of Oxford)-->\n";

	report += "<style type=\"text/css\">\n"\
		  "  table.imagetable {\n"\
		  "  font-family: verdana,arial,sans-serif;\n"\
	  	  "  font-size:12px;\n"\
		  "  color:#333333;\n"\
		  "  border-width: 1px;\n"\
		  "  border-color: #999999;\n"\
		  "  border-collapse: collapse;\n"\
		  "}\n"\
	 	  "table.imagetable th {\n"\
	  	  "  font-size:14px;\n"\
		  "  color: white;\n"
		  "  background:#002147;\n"\
		  "  border-width: 1px;\n"\
		  "  padding: 8px;\n"\
		  "  border-style: solid;\n"\
		  "  border-color: #999999;\n"\
		  "}"\
		  "table.imagetable td {\n"\
		  "  background:#eeeeee;\n"\
		  "  border-width: 1px;\n"\
		  "  padding: 8px;\n"\
		  "  border-style: solid;\n"\
		  "  border-color: #999999;\n"\
		  "}\n"\
		  "</style>\n";

		report += "<table class=\"imagetable\">\n";                          
 	        report += "<tr><th>Equivalent</th>\n";

	if(result.equivalent) {
                report += "<td><b>Yes</b></td></tr>\n";
	} else {
                report += "<td><b>No</b></td></tr>\n";

 	        report += "<tr><th>Counterexample word</th>\n";
    	        for(Word::const_iterator wit=result.ce.begin(); wit!=result.ce.end(); ++wit) {
	          report += "<td> " + alphabet[*wit] + " </td>";
	        }
                report += "  <tr>\n";

		std::string pa = toString(result.a);
		std::string pb = toString(result.b);
		
		report += "<tr><th>Probability in A</th><td>" + pa + "</td></tr>\n";
		report += "<tr><th>Probability in B</th><td>" + pb + "</td></tr>\n";
	}
        report += "</table>";
	report += "<!-- Probabilistic Equivalence Checker (University of Oxford)-->\n";
}

template <typename Automaton>
bool equiv (Automaton& a, Automaton& b, const Settings& settings, Precision<typename Automaton::Weight>& prec) {

  typedef typename Automaton::Matrix Matrix;

  Timer equivalence_timer;

  equivalence_timer.Start();
  
  bool result = false;

  Result<Automaton> r;

  switch(settings.method) {
  case Settings::forward:
	  result = Tzeng<Automaton,true> (a,b,settings,prec,r);
	  break;
  case Settings::backward:
	  result = Tzeng<Automaton,false> (a,b,settings,prec,r);
	  break;
  case Settings::randomised_backward:
	  result = Randomised<Automaton>(a,b,settings,prec,r);
	  break;
  case Settings::randomised_forward:
	   result = Randomised2<Automaton>(a,b,settings,prec,r);
	   break;
  }
  equivalence_timer.Stop();

  

  r.time = equivalence_timer.Read();

  if(result == false && settings.ce) {
    r.a = a(r.ce);
    r.b = b(r.ce);
  }

  if(global_setting::html) {
    std::string report;

    createHTMLReport<Automaton>(r,report);
    std::cout << report << std::endl;
  } else {

	prec.printStatistics();
	if(result == false && settings.ce) {

    
		  std::string word;
		  for(Word::const_iterator wit=r.ce.begin(); wit!=r.ce.end(); ++wit) {
			  word += alphabet[*wit];
			  word += ",";
		  }
		  if(global_setting::verbosity >= 1) {
	 		  std::cout<< "Precision: relative " << prec.relative_precision << " absolute " << prec.absolute_precision << std::endl;
			  std::cout<< "Word \""<<word<<"\" \nA(word) = "<<r.a<<" vs. B(word) = "<< r.b <<" difference " << fabs(r.b - r.a)<<std::endl;
		  }
	  }
    if(global_setting::verbosity >= 1) {
		printf("Time: %5.6f\n",equivalence_timer.Read());
	    printf("Equivalent: %s\n",result ? "yes" : "no");
    }
  }

  return result;
}


#endif
