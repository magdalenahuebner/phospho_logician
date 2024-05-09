% QUERY
% enz_sub_ctxt(Pert, Kpa, Pst, Tprot, Space)

% MAIN RULE (identify context-specific enzyme-substraste relationships)
enz_sub_ctxt(Pert, Kpa, Pst, Tprot, Space):-

  % SEARCH SPACE
  search_space(Space, Kpa, Pst, Tprot), % select 'search space' (pst or tprot; omnipat or pdts or lint)

  % LOGIC MODEL
  enzyme_class(Kpa, EClass, ERes), % source properties: enzyme class
  ksea(Pert, Kpa, Actv, _), % source properties: activity
  perturbs(Pert, Pst, Tprot, Res, Fc, _, _), % target properties: fold change

  check_res(ERes, Res), % check: residue
  check_fc(EClass, Actv, Fc), % check: fold change
  check_sign(Pert, Pst, Tprot, Fc). % check: effect sign


% SEARCH SPACE (either on substrate (pst) level or tprot level)
search_space(omnipath_sub, Kpa, Pst, Tprot):-
  enz_sub(Kpa, Pst, Tprot).

search_space(omnipath_tprot, Kpa, _, Tprot):-
  enz_sub(Kpa, _, Tprot).

search_space(pdts, Kpa, Pst, Tprot):-
  pdt(Kpa, Pst, Tprot).

search_space(pdts_dirts, Kpa, Pst, Tprot):-
  pdt_dirt(Kpa, Pst, Tprot).

search_space(lints, Kpa, Pst, Tprot):-
  us_of(Kpa, Pst, Tprot).

search_space(dirts, Kpa, Pst, Tprot):-
  dirt(Kpa, Pst, Tprot).

search_space(fges, Kpa, Pst, Tprot):-
  fges_enzsub(Kpa, Pst, Tprot).

% RESIDUE (check whether residues targeted by enzyme class and pst residue match)
check_res(ERes, Res):-
  res_match(ERes, Res).

res_match(st, s).
res_match(st, t).
res_match(y, y).
res_match(sty, _).
res_match(unknown, _).


% FOLD CHANGE (check whether source enzyme activity and fold change match)
check_fc(EClass, Actv, Fc):-
  Fc \= unaffected,
  fc_match(EClass, Actv, Fc).

fc_match(kinase, actv, up).
fc_match(kinase, inhb, down).
fc_match(phosphatase, actv, down).
fc_match(phosphatase, inhb, up).


% EFFECT SIGN (check whether target enzyme activity, fc and effect sign match)
check_sign(Pert, Pst, Tprot, _):-
  sign_ignore(Pert, Pst, Tprot).

check_sign(Pert, Pst, Tprot, Fc):-
  ksea(Pert, Tprot, TActv, _), % if target enzyme activity available: assign activity
  effect_sign(Pst, Sign),
  sign_match(TActv, Fc, Sign).

sign_ignore(Pert, Pst, Tprot):-
  \+ ksea(Pert, Tprot, _, _); % if enzyme activity not available (e.g because tprot no kpa)
  ksea(Pert, Tprot, no_t, _); % if enzyme activity not available (e.g because no kpa targets observed)
  ksea(Pert, Tprot, unaffected, _); % if enzyme activity unaffected (other pst might compensate)
  \+ effect_sign(Pst, _); % if effect sign is unavailable
  effect_sign(Pst, unknown); % if effect sign is unknown
  effect_sign(Pst, conflicting). % if effect sign is conflicting

sign_match(actv, up, p_inc).
sign_match(actv, down, p_dec).
sign_match(inhb, up, p_dec).
sign_match(inhb, down, p_inc).
