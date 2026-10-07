load("LMFDB_eigenvalues.sage")

def table_1():
    """
    Compute example in Table 1: mod 3 form with eigenvalues in F9.
    """
    evs = load_all_from_csv(1)
    cands = [h for h in evs if h.p() == 3 and h.level().norm() == 293]
    assert len(cands) == 1
    h = cands[0]
    # assert len(h.get_verified_lifts()) == 1
    h.tabulate_lift_evals()

def table_2():
    """
    Compute example in Table 2: mod 31 form with irrational lift.
    """
    evs = load_all_from_csv(3)
    cands = [h for h in evs if h.p() == 31 and h.level().norm() == 607]
    assert len(cands) == 1
    h = cands[0]
    # assert len(h.get_verified_lifts()) == 1
    h.tabulate_lift_evals()

def table_4():
    """
    Compute example in Table 4: mod 53 form with large char 0 Hecke field.
    """
    evs = load_all_from_csv(7) 
    h = [h for h in evs if len([l for l in h.get_irrational_lifts() if l.coefficient_field().degree() >= 7]) > 0][0]

    return h.tabulate_lift_evals()

def table_5():
    """
    Compute example in Table 5: mod 379 form with cyclotomic Hecke field.
    """
    h = [h for h in load_all_from_csv(7) if h.p() >300 and h.has_lift()][0]
    print(h)
    if not h.has_irrational_lift():
        return "No irrational lifts; rational lifts not currently supported."
    lift = h.get_irrational_lifts()[0]
    F = h.field()
    p = h.p()
    Eself = h.coefficient_field()
    print(f"Table of lifts of {h} over Q(\sqrt(-{h.d()}))")
    E2 = lift.coefficient_field()
    L.<zeta7> = CyclotomicField(7, embedding=None)
    E = L
    i = E2.embeddings(L)[0]    
    cands = [pp for (pp, _) in E.ideal(p).factor() if pp.residue_class_degree() == 1]
    print(f"Splitting behavior of p={p} in Hecke field: (Nm, f, e)",
          [(pp.norm(), pp.residue_class_degree(), pp.ramification_index()) for (pp,_) in E.ideal(p).factor()])
    for pp in cands:
        assert pp.is_principal()
        pp_gen = pp.gens_reduced()[0]
        Epp = E.residue_field(pp)
        for phi in Eself.Hom(Epp):
            good_labels = [lab for lab in lift.good_primes() if ideal_from_label(F, lab).is_coprime(F.ideal(p))]
            if all([phi(e) == Epp(i(lift.eigenvalues()[lab])) for (lab,e) in h.eigenvalues().items() if lab in good_labels]):
                print(f"Found correct ideal, reducing mod {pp_gen}")
                print(f"Lift has level {lift.level_label()} = {[(ideal_label(pp), m) for (pp,m) in lift.level().factor()]}")
                print("Hecke eigenvalues defined over", E2, "of discriminant", E2.discriminant(), "=", E2.discriminant().factor())
                table = [[lab + ("*" if lab not in good_labels else ""),
                            h.eigenvalues()[lab],
                            Epp(i(lift.eigenvalues()[lab])),
                            E(i(lift.eigenvalues()[lab])),
                            ] for lab in h.eigenvalues()]
                print(tabulate(table, tablefmt="simple_grid", headers=["$\p$",
                                                                  "$a_{\p}(f)$",
                                                                  f"$a_{{\p}}(F) \Mod {pp_gen}$",
                                                                  "$a_{\p}(F)$",
                                                                  ]))
    return 1

    

def non_steinberg_lift():
    """
    Compute Example 4.8: non-Steinberg lift at p=2.
    """
    from lmf import db
    cands = [h for h in load_liftable_from_csv(11) if h.p() == 2 and h.level_label() == "23.1"]
    assert len(cands) == 1
    h = cands[0]
    print("Number of lifts:", len(h.get_rational_lifts()))
    lift_reduction_types(h)

def table_6():
    """
    Compute example in Table 6: Multiplicity 2 form of level 931.2
    """
    print("This form is not stored because the eigenvalue system lifts. Compute directly in Magma instead.")
    
def table_7(num_primes=11, good_primes_only = True):
    """
    Compute example in Table 7: dihedral lift along with splitting behaviour in extension.
    """
    evs = load_all_from_csv(7)
    cands = [h for h in evs if h.p() == 7 and h.level().norm() == 512]
    assert len(cands) == 1
    h = cands[0]
    K = h.field()
    w = K.gen()
    assert w^2 - w + 2 == 0
    R.<x> = PolynomialRing(K)
    L.<z> = K.extension(x^2+w)
    if good_primes_only:
        keys = h.good_primes()[:num_primes]
    else:
        keys = list(h.eigenvalues())[:num_primes]

    # Helper function printing splitting of prime ideal in quadratic extension L
    def quad_splt(pp):
        fac = L.ideal(pp).factor()
        if len(fac) == 1:
            if fac[0][1] > 1:
                return "ram"
            return "inert"
        return "split"

    table = []
    table.append(["pp"] + keys)
    table.append(["pp*OL"] + [quad_splt(ideal_from_label(K, key)) for key in keys])
    table.append(["a_pp(f)"] + [h.eigenvalues()[key] for key in keys])
    print(tabulate(table, tablefmt = "simple_grid"))
    
    
def table_8():
    """
    Compute table of Fp-valued systems, admissible and not, with LMFDB lifts
    """
    table = []
    for d in [1,2,3,7,11]:
        evs = load_all_from_csv(d)
        all_p = [h for h in evs if h.is_rational()]
        small_p = [h for h in all_p if h.p() < 20]
        assert len([h for h in small_p if not h.is_admissible() and h.has_rational_lift()]) == 0
        table.append([
            d,
            len(small_p),
            len([h for h in small_p if h.is_admissible()]),
            len([h for h in small_p if h.has_rational_lift()]),
            len(all_p),
            len([h for h in all_p if h.is_admissible()]),
            len([h for h in all_p if h.has_rational_lift()]),
        ])
    print(tabulate(table, tablefmt="simple_grid", headers=["d", "Eth p<20", "Adm", "LMFDB ✓",
                                                           "Eth all p", " Adm", "LMFDB ✓"]))


def table_9():
    """
    Compute table of lifts of ethereal BMFs (not restricted to LMFDB)
    """
    all_p_rows = []
    some_p_rows = []
    for d in [1,2,3,7,11]:
        all_p_forms = load_all_from_csv(d)
        some_p_forms = [h for h in all_p_forms if h.p() < 20]
        all_p_rows.append([
            d,
            len(all_p_forms),
            len([h for h in all_p_forms if h.has_lift()]),
            len([h for h in all_p_forms if h.get_verified_lifts()]),
        ])
        some_p_rows.append([
            d,
            len(some_p_forms),
            len([h for h in some_p_forms if h.has_lift()]),
            len([h for h in some_p_forms if h.get_verified_lifts()]),
        ])
    print("Table of ethereal forms with p < 20")
    print(tabulate(some_p_rows,tablefmt="simple_grid", headers=["d", "Eth", "Lifts", "Verif"]))
    print("Table of ethereal forms, all p")
    print(tabulate(all_p_rows, tablefmt="simple_grid", headers=["d", "Eth", "Lifts", "Verif"]))

#############################
# Helper functions
#############################

def lift_reduction_types(h, only_new_level=True, kodaira_all_curves=False):
    """Print table of rational lifts along with reduction types at bad primes.

    The reduction type and conductor exponent only depend on the isogeny class.
    The Kodaira symbol depends on the individual curve, so by default we print
    the one corresponding to the first curve in the class, per LMFDB labeling.
    Set kodaira_all_curves to 'True' to list it for all of them.
    """
    reduction_dict = {1: 'split mult (Steinberg)',
           -1: 'nonsplit mult (Steinberg)',
           0: 'additive (non-Steinberg)',
           }

    from lmf import db
    K = h.field()
    labels = h.get_rational_lifts()
    query = {'class_label': {'$in': labels}}
    if not kodaira_all_curves:
        query['number'] = 1
    found = {}
    for r in db.ec_nfcurves.search(query,
            ['class_label', 'conductor_label', 'local_data', 'semistable', 'number']):
        found.setdefault(r['class_label'], []).append(r)

    table = []
    for lab in labels:
        curves = found.get(lab)
        if curves is None:
            table.append([lab, '-', 'no curve in LMFDB', '', '', ''])
            continue
        kods = {}
        for c in curves:
            for ld in c['local_data']:
                q = K.ideal(sage_eval(ld['p'], locals={'w': K.gen()}))
                kods.setdefault(q, []).append(str(KodairaSymbol(ld['kod'])))
        for ld in curves[0]['local_data']:
            q = K.ideal(sage_eval(ld['p'], locals={'w': K.gen()}))
            if only_new_level and q.divides(h.level()):
                continue                      # skip primes already in the mod p level
            table.append([lab, prime_label(q), reduction_dict[ld['red']],
                          ", ".join(dict.fromkeys(kods[q])),
                          ld['ord_cond'], curves[0]['semistable']])
    print(tabulate(table,
                   headers=['lift', 'q', 'reduction at q', 'Kodaira', 'ord_q(N)', 'semistab'],
                   tablefmt='simple_grid'))
    return table

def theorem_counts():
    """
    Compute number of liftable forms in total, as in Theorem 1.3 and Theorem 4.3.
    """
    evs = load_forms_all_d()
    all_lifties = [h for h in evs if h.has_lift()]
    print(f"Number of liftable forms: {len(all_lifties)}/{len(evs)}= {n(len(all_lifties)/len(evs))}.")
    restr_evs = [h for h in evs if 2 < h.p() < 20 and h.level().norm() <= 500]
    restr_lifties = [h for h in restr_evs if h.has_lift()]
    print(f"Number of liftable forms Nm(n) < 500, 2 < p < 20: {len(restr_lifties)}/{len(restr_evs)} = {n(len(restr_lifties)/len(restr_evs))}")
    
