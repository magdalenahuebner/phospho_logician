direct_target(K, Pst, Tprot):-
  us_of(K, Pst, Tprot),
  findall(X, (
    us_of(X, _, Tprot), 
    X \= K, 
    X \= Tprot,
    us_of(K, _, X),
    \+ us_of(Tprot, _, X)
    ), List),
  List = [].

direct_target(K, Pst, Tprot, no_loop):-
  us_of(K, Pst, Tprot),
  findall(X, (
  	us_of(X, _, Tprot), 
  	X \= K, 
  	X \= Tprot, 
  	us_of(K, _, X)
  	), List),
  List = [].

direct_target(K, Pst, Tprot, loop):-
  us_of(K, Pst, Tprot),
  findall(X, (
    us_of(X, _, Tprot), 
  	X \= K, 
  	X \= Tprot, 
  	us_of(K, _, X),
  	us_of(Tprot, _, X)
  	), List),
  List \= [].

indirect_target(Pert, K, Pst):-
  knowninhibitor(Pert, K),
  perturbs(Pert, Pst, _, Fc, _, _),
  Fc \= unaffected.



/*
% test
us_of(c, 1, f).
us_of(c, 2, f).
us_of(c, 1, g).
us_of(c, 1, k).
us_of(c, 1, l).
us_of(c, 1, m).
us_of(c, 1, o).
us_of(c, 1, p).
us_of(c, 1, q).
us_of(f, 1, k).
us_of(f, 1, o).
us_of(f, 1, l).
us_of(f, 1, p).
us_of(k, 1, l).
us_of(k, 1, p).
us_of(k, 1, o).
us_of(l, 1, p).
us_of(g, 1, m).
us_of(g, 1, q).
us_of(g, 1, l).
us_of(g, 1, p).
us_of(m, 1, q).
us_of(m, 1, l).
us_of(m, 1, p).
us_of(m, 1, g).
us_of(d, 1, g).
us_of(d, 2, g).
us_of(d, 1, m).
us_of(d, 1, q).
us_of(d, 1, l).
us_of(d, 1, p).
*/

