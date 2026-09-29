# Uses lmfdb lite, see https://github.com/roed314/lmfdb-lite
# Make sure to install via sage's own pip:
# sage -pip install -U "lmfdb-lite[pgbinary] @ git+https://github.com/roed314/lmfdb-lite.git@4c4262b60b6ae8a9c6f33771ae4e939637d95741"
# and tabulate:
# sage -pip install -U tabulate

from tabulate import tabulate
# Stolen from lmfdb - used for sorting ideals
from psort import primes_iter, prime_label, prime_from_label, ideal_label, ideal_from_label
import re
import csv

SEARCH_COUNT = ZZ(1e6)

# NOTE: We cannot search on Hecke eigenvalue modulo p because the 'hecke_eigs' field is jsonb, not array.
# So we have to search through all the forms. 

def QuadFld(d):
    """
    Helper function to create quadratic number
    field with monogenic integer ring}
    """
    if (d-3) % 4 == 0:
        pol = x^2-x+(d+1)/4
    else:
        pol = x^2+d

    F.<t> = NumberField(pol)
    ZF = F.maximal_order()

    return F

class HeckeEig():
    def __init__(self, evals, lvl, p, label=None, eig_dim = None, gen_dim = None, source = None):
        self.evals = evals
        self.lvl = lvl
        self.prime = p
        self.rational_lifts = []
        self.irrational_lifts = []
        self._label = label
        self.eig_dim = eig_dim
        self.gen_dim = gen_dim
        # For a char. 0 lift: the (level, p) of the mod p form it lifts
        # in ../data/irrat_lifts_d*.csv.
        self.source = source

    def conjugate_form(self):
        """
        Initialize conjugate form (WARN: does not copy lifts or other invariants)
        """
        K = self.field()
        c = [c for c in K.automorphisms() if c(K.gen()) != K.gen()][0]
        conj_label = lambda label: ideal_label(c(ideal_from_label(K, label)))
        return HeckeEig(
            {lab: self.eigenvalues()[conj_label(lab)] for lab in self.eigenvalues()},
            c(self.level()),
            self.p(),
        )

    def __repr__(self):
        if self._label:
            return self.label()
        else:
            return f"Characteristic {self.p()} Bianchi eigenvalue system of level {self.level_label()}"

    # Helper methods to make code below cleaner.
    # Not strictly necessary, but prevents us from accidentally
    # changing variables in HeckeEig object.
    def level(self):
        return self.lvl

    def level_label(self):
        return ideal_label(self.level())
    def label(self):
        return self.label
    
    def field(self):
        return self.level().number_field()

    def p(self):
        return self.prime

    def d(self):
        return squarefree_part(-self.field().discriminant())
    
    def eigenvalues(self):
        return self.evals

    def eigenvalue_list(self):
        return [v for _,v in self.eigenvalues().items()]

    def add_rational_lift(self, lift):
        self.rational_lifts.append(lift)

    def add_irrational_lift(self, lift):
        self.irrational_lifts.append(lift)

    def has_rational_lift(self):
        return self.rational_lifts != []

    def has_irrational_lift(self):
        return self.irrational_lifts != []

    def has_lift(self):
        return self.has_rational_lift() or self.has_irrational_lift()

    def get_rational_lifts(self):
        return self.rational_lifts

    def get_verified_lifts(self):
        lift_label = lambda lift: ideal_label(level_from_BMF_label(self.field(), lift))
        # Since verified_d*.csv uses old data, need to check conjugate level of lift as well.
        c = self.field().automorphisms()[1]
        get_conj_label = lambda label: ideal_label(c(ideal_from_label(self.field(), label)))

        conj_form = self.conjugate_form()
        p_str = str(self.p())
        
        evs_match = lambda evs1, evs2: all((ev == evs2[i] or evs2[i] == 0 for i,ev in enumerate(evs1) ))

        verified_lifts = []
        with open(f"../data/verified_d{self.d()}.csv") as f:
            csv_reader = csv.reader(f, delimiter=';')
            for l in csv_reader:
                file_evs = eval(l[3])
                cand = self if l[0] == self.level_label() else conj_form;
                if l[0] == cand.level_label() and l[1] == p_str and evs_match(cand.eigenvalue_list(), file_evs):
                    lift = next((lift for lift in self.get_rational_lifts() if
                                 l[2] == (lift_label(lift) if cand is self else get_conj_label(lift_label(lift)))), None)
                    if lift is not None:
                        verified_lifts.append(lift)
        return verified_lifts

    def get_irrational_lifts(self):
        return self.irrational_lifts

    def get_lifts(self):
        return self.get_rational_lifts() + self.get_irrational_lifts()

    def source_key(self):
        return self.source

    def is_rational(self):
        return self.eigenvalue_list()[0].parent().degree() == 1

    def eigenspace_dim(self):
        if self.eig_dim:
            return self.eig_dim
        else:
            return -1

    def coefficient_field(self):
        return self.eigenvalue_list()[0].parent()

    def generalized_eigenspace_dim(self):
        if self.gen_dim:
            return self.gen_dim
        else:
            return -1

    def good_primes(self):
        p = self.p()
        return [lab for lab in self.eigenvalues() if ideal_from_label(self.field(), lab).is_coprime(self.level() * (p if p != 0 else 1))]

    def tabulate_lift_evals(self, tablefmt="simple_grid", horizontal=False):
        """
        Return a table with the eigenvalues of form and its lifts
        """
        if self.get_irrational_lifts() == []:
            return "No irrational lifts; rational lifts not currently supported."
        lifts = self.get_lifts()
        F = self.field()
        p = self.p()
        Eself = self.coefficient_field()
        print(f"Table of lifts of {self} over Q(\sqrt(-{self.d()}))")
        for lift in lifts:
            E = lift.coefficient_field()
            # for pretty printing:
            E.<w> = E.change_names()
            cands = [pp for (pp, _) in E.ideal(p).factor() ]
            print(f"Splitting behavior of p={p} in Hecke field: (Nm, f, e)", [(pp.norm(), pp.residue_class_degree(), pp.ramification_index()) for (pp,_) in E.ideal(p).factor()])
            for pp in cands:
                assert pp.is_principal()
                pp_gen = pp.gens_reduced()[0]
                Epp = E.residue_field(pp)
                for phi in Eself.Hom(Epp):
                    good_labels = [lab for lab in lift.good_primes() if ideal_from_label(F, lab).is_coprime(F.ideal(p))]
                    if all([phi(e) == Epp(lift.eigenvalues()[lab]) for (lab,e) in self.eigenvalues().items() if lab in good_labels]):
                        print(f"Found correct ideal, reducing mod {pp_gen}")
                        print(f"Lift has level {lift.level_label()} = {[(ideal_label(pp), m) for (pp,m) in lift.level().factor()]}")
                        print("Hecke eigenvalues defined over", E, "of discriminant", E.discriminant(), "=", E.discriminant().factor())
                        table = [[lab + ("*" if lab not in good_labels else ""),
                                  self.eigenvalues()[lab],
                                  Epp(lift.eigenvalues()[lab]),
                                  E(lift.eigenvalues()[lab]),
                                  ] for lab in self.eigenvalues()]
                        print(tabulate(table, tablefmt=tablefmt, headers=["$\p$",
                                                                   "$a_{\p}(f)$",
                                                                   f"$a_{{\p}}(F) \Mod {pp_gen}$",
                                                                   "$a_{\p}(F)$",
                                                                          ]))

    def is_CM(self):
        """Check if form is (likely CM) by checking if enough eigenvalues are 0."""
        evals = list(self.eigenvalues().values())
        return len([x for x in evals if x == 0])/len(evals) > 1/3
        
    def csv_dict(self):
        if self.has_rational_lift():
            return {'level': self.level_label(),
                    'evals': self.eigenvalues(),
                    'p': self.p(),
                    'eig_dim': self.eigenspace_dim(),
                    'gen_dim': self.generalized_eigenspace_dim(),
                    'lift_labels': self.get_rational_lifts(),
                    }
        else:
            return {'level': self.level_label(),
                    'evals': self.eigenvalues(),
                    'p': self.p(),
                    'eig_dim': self.eigenspace_dim(),
                    'gen_dim': self.generalized_eigenspace_dim(),
                    'Coeff_minpoly': self.eigenvalue_list()[0].parent().gen().minpoly(),
                    }

    def non_LR_primes(self, lift):
        """
        return list of primes not satisfying the level raising condition
        """
        non_LR_primes = []
        p = self.p()
        prime_list = primes_iter(self.field())
        for i, pp in enumerate(prime_list):
            if i >= len(self.eigenvalues()):
                return non_LR_primes
            elif pp.divides(lift['level']) and not pp.divides(p*self.level()):
                modp_diff = (self.eigenvalues()[i]^2 - (1 + norm(pp))^2) % p
                if modp_diff != 0:
                    non_LR_primes.append(prime_label(pp))

    def non_LR_lifts(self):
        """
        return dict of pairs lift_label: [ideal_labels] which fail the level raising condition
        """
        lift_dict = {}
        for lift in self.get_rational_lifts():
            lift_dict[lift['label']] = self.non_LR_primes(lift)
        return lift_dict
    
    def LR_primes(self, bd=100):
        LR_primes = []
        p = self.p()
        for i,pp in enumerate(primes_iter(self.field())):
            if pp.norm() > bd or i >= len(self.eigenvalues()):
                return LR_primes
            else:
                modp_diff = (self.eigenvalues()[i]^2 - (1 + norm(pp))^2) % p
                if modp_diff == 0:
                    LR_primes.append(prime_label(pp))
        
        

def find_congruent_forms(d, input_file = None, find_all_lifts = False):
    """
    Search LMFDB for Bianchi eigenforms on field on discriminant d
    with eigenvalues congruent to eigenvalue tuples stored in `../data/nonEis_d*.csv`.
    Returns list of hecke eigenvalue objects, whose field `rational_lifts` is populated
    with rational lifts, in the form of a tuple (label, evs, level)
    """
    # Place import here to improve loading speed of script
    from lmf import db

    R.<x> = PolynomialRing(ZZ)
    K = QuadFld(d)
    BMFs = db.bmf_forms
    OK = K.maximal_order()
    label = f"2.0.{K.discriminant().abs()}.1"
    if input_file is None:
        input_file = f"../data/nonEis_d{d}.csv"

    # load mod p BMFs from input_file
    # this is small enough that it doesn't make sense to use an iterator.
    print(f"Loading eigenvalue systems from {input_file}, this might take some time.")
    # modp_evs = [h for h in load_nonEis(d) if h.is_rational()]
    modp_evs = load_nonEis(d)
    max_ev_len = max(len(h.eigenvalues()) for h in modp_evs)
    # extract keys from eigenvalue system
    prime_list = [prime_from_label(K,lab) for lab in modp_evs[0].eigenvalues()]
    print([ideal_label(pp) for pp in prime_list])
    print(f"Successfully loaded {len(modp_evs)} eigenvalues from {input_file}")

    hits = BMFs.search({'dimension': 1, "field_label": label}, ['label', 'hecke_eigs', 'level_bad_primes', 'level_label'], limit = SEARCH_COUNT)
    print("Successfully pulled Hecke eigenvalue systems from LMFDB")

    for hit in hits:
        hit_lvl = ideal_from_label(K, hit['level_label'])
        hit_ev = hit['hecke_eigs']

        for h in modp_evs:
            # Irrational forms can't have rational lifts;
            # unless we're trying to find all lifts, skip to next once we have a rational lift
            if not h.is_rational() or (not find_all_lifts and h.has_rational_lift()):
                continue

            lvl = h.level()

            # if the level of potential lift isn't divisible by level of hecke_ev, skip to next
            if not lvl.divides(hit_lvl):
                continue
            # need list of eigenvalues to match up with LMFDB results
            ev = list(h.eigenvalues().values())
            p = h.p()

            min_len = min(len(ev), len(hit_ev))
            # We allow a lift to have different eigenvalues corresponding to "bad" primes
            # This includes:
            #  - Primes dividing p,
            #  - Primes above the level of the lift
            #  - Primes dividing the discriminant of the field

            if all(((ev[i] - hit_ev[i]) % p == 0 or prime_list[i].divides(p*hit_lvl))
                   for i in range(min_len)):
                print("Found rational lift for", h, "of level", ideal_label(hit_lvl))
                h.add_rational_lift(hit['label'])
            
    return modp_evs

def test(d, find_all_lifts=False):
    modp_evs = find_congruent_forms(d, find_all_lifts=find_all_lifts)

    for hecke_ev in modp_evs:
        if hecke_ev.has_rational_lift():
            evs = hecke_ev.eigenvalues()
            print(f"Mod {hecke_ev.p()} level {hecke_ev.level_label()} eigenvalue system {hecke_ev.eigenvalues()} has lifts:")
            for lift in hecke_ev.get_rational_lifts():
                print(lift['label'], "--- level:", ideal_label(lift['level']) )  # label
                print(lift['evals'][:len(evs)]) # first few eigenvalues
                p = hecke_ev.p()
                print("Difference between eigenvalues mod p:", [(lift['evals'][i] - evs[i]) % p for i in range(len(evs))])
            print("-"*10, "\n")

    print("Hecke eigenvalue systems with no lifts:")
    print(*(x for x in modp_evs if not x.has_rational_lift()), sep='\n')
    return modp_evs

def write_congruent_forms(d, modp_evs = None):
    if modp_evs is None: 
        modp_evs = find_congruent_forms(d, find_all_lifts=true)

    with open(f'd{d}_liftable.csv', 'w', newline='') as f:
        fieldnames = ['level','p','evals','lift_labels','eig_dim','gen_dim']
        writer = csv.DictWriter(f, delimiter=';', fieldnames=fieldnames)
        writer.writeheader()
        for h in modp_evs:
            if h.has_rational_lift():
                writer.writerow(h.csv_dict())
    print(f"Wrote liftable eigenvalues and lifts to", f'd{d}_liftable.csv')

    with open(f'd{d}_unliftable.csv', 'w', newline='') as f:
        fieldnames = ['level','p','evals','eig_dim','gen_dim', 'Coeff_minpoly']
        writer = csv.DictWriter(f, delimiter=';', fieldnames=fieldnames)
        writer.writeheader()
        for h in modp_evs:
            if not h.has_rational_lift():
                writer.writerow(h.csv_dict())
    print(f"Wrote unliftable eigenvalues to", f'd{d}_unliftable.csv')

def load_nonEis(d):
    """ Load all Hecke eigenvalue systems from nonEis_dX.csv files
    """
    K = QuadFld(d)
    
    input_file = f"../data/nonEis_d{d}.csv"
    
    modp_evs = []
    with open(input_file) as f:
        csv_reader = csv.reader(f, delimiter=';')

        # iterate once to get past header line(?)
        # print(next(csv_reader))
        prime_labels = next(csv_reader)[2:-3]

        _ = magma.eval('R<x> := PolynomialRing(Integers())')

        for row in csv_reader:
            print(row)
            if row == [] or "ERROR" in row[2]:
                continue
            lvl = ideal_from_label(K, row[0])
            p = ZZ(row[1])
            min_poly = row[-1]
            gen_dim = ZZ(row[-2])
            eig_dim = ZZ(row[-3])
            evs = {}
            is_valid = True
            F = magma(f"ext<GaloisField({p}) | {min_poly}>")
            for i, coeff in enumerate(row[2:-3]):
                try:
                    evs[prime_labels[i]] = F(coeff).sage()
                    
                except TypeError:
                    if coeff == "-":
                        evs[prime_labels[i]] = (ZZ(0))
                    else:
                        print(f"Failed to parse {row} ")
                        is_valid = False
                        break
            if is_valid:
                modp_evs.append(HeckeEig(
                    evs, lvl, p,
                    eig_dim=eig_dim,
                    gen_dim=gen_dim,
                ))
    return modp_evs
    
def load_unliftable_from_csv(d):
    K = QuadFld(d)
    modp_evs = []
    with open(f'd{d}_unliftable.csv', 'r', newline='') as f:
        # fieldnames = ['level','p','evals','lift_labels', 'lift_evs']
        reader = csv.DictReader(f, delimiter=';')
        for r in reader:
            lvl = ideal_from_label(K, r['level'])
            p = ZZ(r['p'])
            ZZx.<x> = PolynomialRing(ZZ)
            E.<w> = FiniteField(p).extension(sage_eval(r['Coeff_minpoly'], locals={'x':x}))
            eigs = sage_eval(r['evals'], locals={'X1':w})
            eigs = {key: E(e) for (key, e) in eigs.items()}

            modp_evs.append(HeckeEig(eigs, lvl, p, eig_dim=ZZ(r['eig_dim']), gen_dim=ZZ(r['gen_dim'])))
    return modp_evs

def load_liftable_from_csv(d):
    K = QuadFld(d)
    modp_evs = []
    with open(f'd{d}_liftable.csv', 'r', newline='') as f:
        # fieldnames = ['level','p','evals','lift_labels', 'lift_evs']
        reader = csv.DictReader(f, delimiter=';')
        for r in reader:
            p = ZZ(r['p'])
            lvl = ideal_from_label(K, r['level'])
            E = FiniteField(p)
            eigs = sage_eval(r['evals'])
            eigs = {key: E(e) for (key, e) in eigs.items()}
            h = HeckeEig(eigs, lvl, p, eig_dim=ZZ(r['eig_dim']), gen_dim=ZZ(r['gen_dim']))
            for lift in eval(r['lift_labels']):
                h.add_rational_lift(lift)
            modp_evs.append(h)
            
    return modp_evs

def load_irrat_lifts(d, modp_evs=None):
    """Load list of all irrational forms with lifts from ../data/irrat_lifts_d{d}.csv,
    ignoring the ones where we haven't found anything yet.
    If modp_evs is not None but a list of eigenvalues, then record the lifts.

    Since this is slow, we cache the results.
    """
    K = QuadFld(d)
    irrat_lifts = []
    irrat_file= f'../data/irrat_lifts_d{d}.csv'
    print(f"Loading irrational lifts from {irrat_file}, this might take a while.")
    _ = magma.eval('R<x> := PolynomialRing(Integers())')

    # extract prime labels first, since to access them in DictReader
    with open(irrat_file) as f:
        reader = csv.reader(f, delimiter=';')
        prime_labels = next(reader)[3:-1]
    
    with open(irrat_file) as f:
        csv_reader = csv.DictReader(f, delimiter=';')
        for r in csv_reader:
            if not "No" in r[" lift level"]:
                E = magma(f'NumberField({r["Coeff_minpoly"]})')
                lvl = ideal_from_label(K, r["level"])
                lift_lvl = ideal_from_label(K, r[" lift level"])
                p = ZZ(r[" prime"])
                evals = {lab: E(r[lab]).sage() for lab in prime_labels}
                h = HeckeEig(evals, lift_lvl, 0, source=(lvl, p))
                irrat_lifts.append(h)
    if not modp_evs:
        modp_evs = load_unliftable_from_csv(d)

    lifts_by_key = {}
    for h in irrat_lifts:
        lifts_by_key.setdefault(h.source_key(), []).append(h)

    forms_by_key = {}
    for g in modp_evs:
        forms_by_key.setdefault((g.level(), g.p()), []).append(g)

    undetected = []
    for key, lifts in lifts_by_key.items():
        forms = forms_by_key.get(key, [])
        if len(forms) == 1:
            # Unique lift coming from this level, so can add directly
            for h in lifts:
                forms[0].add_irrational_lift(h)
            continue
        for g in forms:
            # ...otherwise, match up eigenvalues
            
            match = next((h for h in lifts if is_irrational_lift(h, g)), None)
            if match is None:
                continue
            g.add_irrational_lift(match)
            lifts.remove(match)
        undetected += lifts
    if len(undetected) > 0:
        print("WARNING: did not find mod p forms for", undetected)
    return irrat_lifts

def is_irrational_lift(h, g):
    """
    Return True if h is an irrational char. 0 lift of g.
    """
    p = g.p()
    Eg = g.coefficient_field()
    E = h.coefficient_field()
    for pp,_ in factor(E.ideal(p)):
        Emodp = pp.residue_field()
        if Emodp.cardinality() <= Eg.cardinality():
            for phi in Eg.Hom(Emodp):
                good_primes = [pp for pp in h.good_primes() if pp in g.good_primes()]
                matches = [phi(g.eigenvalues()[lab]) == Emodp(h.eigenvalues()[lab]) for lab in good_primes]
                if all(matches):
                    return True
    return False
            


def load_all_from_csv(d):
    """
    Load liftable and unliftable forms over field Q(sqrt{-d}), including those where we haven't found lifts yet.
    """
    lifties = load_liftable_from_csv(d)
    unlifties = load_unliftable_from_csv(d)
    _ = load_irrat_lifts(d, modp_evs = unlifties)
    
    return lifties + unlifties

def load_forms_all_d():
    """
    Load all forms across all fields.
    """
    return sum([load_all_from_csv(d) for d in [1,2,3,7,11]], [])
    
def ul_filter_csv(d):
    """
    Create magma-readable file with EBMFs which don't lift in the lmfdb
    """
    F = QuadFld(d)
    input_file = f"../data/nonEis_d{d}.csv"
    output_file = f"../data/lmfdbNonlift_d{d}.csv"

    # irrational forms automatically don't lift, so copy line from nonEis directly
    modp_evs = [h for h in load_unliftable_from_csv(d) if h.coefficient_field().cardinality().is_prime()]
    lines = []
    with open(input_file, newline='') as f:
        csv_reader = csv.reader(f, delimiter=';')
        # add top line
        lines.append(";".join(next(csv_reader)))
        for i,row in enumerate(csv_reader):
            if i == []:
                continue
            p = eval(row[1])
            # irrational eigenvalue systems are never searched for in the LMFDB,
            # so they always get passed on to the magma stage
            if any("[" in r for r in row):
                lines.append(";".join(row))
                continue
            evs = []
            is_valid = True
            for x in row[2:-3]:
                evs.append(ZZ(x))
            cands = [h for h in modp_evs if all([x == evs[i]  for (i,x) in enumerate(h.eigenvalue_list())])]
            if cands:
                lines.append(";".join(row))
                
            # if any(row[0] == h.level_label() and p == h.p() and h.eigenvalue_list() == evs for h in modp_evs):

    
    with open(output_file, "w") as f:
        f.write("\n".join(lines))
    print(f"Wrote unliftable forms to {output_file}")

def ul_print_LR_primes(d):
    F = QuadFld(d)
    # modp_evs = [ev for ev in find_congruent_forms(d) if not ev.has_rational_lift()]
    modp_evs = load_unliftable_from_csv(d)
    modp_evs = [h for h in modp_evs if h.LR_primes()]
    
    for h in modp_evs:
        print(h.level_label(),h.p(),h.LR_primes())
        
def print_lifts(d):
    modp_evs = load_liftable_from_csv(d)
    F = QuadFld(d)
    for h in modp_evs:
        # for lab in h.get_rational_lifts():
        #     print(lab)
        lift_levels = [level_from_BMF_label(F, lab) for lab in h.get_rational_lifts()]
        print(f"Mod {h.p()} eigenform of lvl {h.level_label()} has lifts from raising level at ideals", *[ideal_label(lvl/h.level()) for lvl in lift_levels])
        
def level_from_BMF_label(F, label):
    return ideal_from_label(F, label.split('-')[1])

def load_p_liftable(d):
    evs = load_liftable_from_csv(d)
    F = QuadFld(d)
    p_evs = []
    p_count = 0
    for ev in evs:
        has_p_lift = False
        p = ev.p()
        lvl = ev.level()
        for lift in ev.get_rational_lifts():
            lift_lvl = level_from_BMF_label(F, lift)
            factors = [fac[0] for fac in (lift_lvl/lvl).factor()]
            if not any([fac.is_coprime(p) for fac in factors]):
                print(p, ideal_label(lvl), ideal_label(lift_lvl/lvl))
                p_evs.append(ev)
                break
    print(f"{len(p_evs)}/{len(evs)} have p-lifts")
    return p_evs
        
def count_l_vs_ul(all_evs, norm_bd=500, p_bd=20, print_unliftable=False):
    ls = [h for h in all_evs if h.has_lift()]
    uls = [h for h in all_evs if not h.has_lift()]
    ls = [ev for ev in ls if ev.level().norm() <= norm_bd and ev.p() < p_bd]
    uls = [ev for ev in uls if ev.level().norm() <= norm_bd and ev.p() < p_bd]
    print(f"{len(ls)}/{len(ls) + len(uls)} of forms with d = {ls[0].field().discriminant()}, p < {p_bd} and lvl norm < {norm_bd} have char 0 lifts")
    if print_unliftable and len(uls) > 0:
        print("Forms without known lift:" + " "*(len(f"{uls[0]}") - 18) + "Degree of Hecke field")
        for ev in uls:
            print(f"{ev}\t {ev.coefficient_field().degree()}")
    return len(ls)/(len(ls) + len(uls))

    
def count_single_prime_lifts(d):
    F = QuadFld(d)
    ls = load_liftable_from_csv(d)
    def is_good(h):
        newlvl = level_from_BMF_label(F, h.get_rational_lifts()[0])/h.level()
        return newlvl.is_prime() and newlvl.is_coprime(h.p())
    return len([h for h in ls if is_good(h)])/len(ls)

def liftable_statistics():
    def make_table(forms):
        return [["Ethereal forms total:", len(forms)],
                    ["Ethereal forms with lift:", len([h for h in forms if h.has_lift()])],
                    ["Ethereal Forms with no lift:", len([h for h in forms if not h.has_lift()])],
                    ["Fp-valued forms total:", len([h for h in forms if h.is_rational()])],
                    ["Fp-valued forms with rational lift:", len([h for h in forms if h.is_rational() and h.has_rational_lift()])],
                    ["Fp-valued forms with irrational lift:", len([h for h in forms if h.is_rational() and h.has_irrational_lift()])],
                    ["Fp-valued forms with lift (total):", len([h for h in forms if h.is_rational() and h.has_lift()])],
                    ["Fp-valued forms with no lift:", len([h for h in forms if not h.has_lift() and h.is_rational()])],
                ]
    total_forms = []
    total_small_forms = []
    total_small_p_forms = []
    for d in [1,2,3,7,11]:
        all_forms = load_all_from_csv(d)
        small_p_forms = [h for h in all_forms if h.p() < 20]
        small_forms = [h for h in all_forms if h.p() < 20 and h.level().norm() <= 500]
        total_forms += all_forms
        total_small_p_forms += small_p_forms
        total_small_forms += small_forms
        print("\n")
        print(tabulate(make_table(all_forms), headers=[f"d={d}, no conditions", "Count"], tablefmt="grid"), "\n")
        print(tabulate(make_table(small_p_forms), headers=[f"d={d}, p < 20", "Count"], tablefmt="simple_grid"))
        print(tabulate(make_table(small_forms), headers=[f"d={d}, p < 20, Nm(n) <= 500", "Count"], tablefmt="simple_grid"))
        print("\n")
    print("\n")
    print(tabulate(make_table(total_forms), headers=[f"All d, no conditions", "Count"], tablefmt="double_grid"), "\n")
    print(tabulate(make_table(total_small_p_forms), headers=[f"All d, p < 20", "Count"], tablefmt="double_grid"))
    print(tabulate(make_table(total_small_forms), headers=[f"All d, p < 20, Nm(n) <= 500", "Count"], tablefmt="double_grid"))
    print("\n")
    eigen_table = [["Ethereal forms with multiplicity >1:", len([h for h in all_forms if h.eigenspace_dim() > 1])],
                   ["Ethereal forms with extra generalized eigenspace:", len([h for h in all_forms if h.generalized_eigenspace_dim() > h.eigenspace_dim()])],
                   ["Ethereal forms with extra generalized eigenspace, p>3:", len([h for h in all_forms if h.generalized_eigenspace_dim() > h.eigenspace_dim() and h.p() > 3])],
                   ]
    print(tabulate(eigen_table))

def table_9():
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
    

def big_prime_example(tablefmt="simple_grid"):
    """
    Print eigenvalue table for mod 379 lift, with char 0 eigenvalues coerced to cyclotomic field
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
    print(f"Splitting behavior of p={p} in Hecke field: (Nm, f, e)", [(pp.norm(), pp.residue_class_degree(), pp.ramification_index()) for (pp,_) in E.ideal(p).factor()])
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
                print(tabulate(table, tablefmt=tablefmt, headers=["$\p$",
                                                                  "$a_{\p}(f)$",
                                                                  f"$a_{{\p}}(F) \Mod {pp_gen}$",
                                                                  "$a_{\p}(F)$",
                                                                  ]))
    return 1


def big_field_example():
    """
    Print eigenvalue table of large degree
    """
    evs = load_all_from_csv(7) 
    h = [h for h in evs if len([l for l in h.get_irrational_lifts() if l.coefficient_field().degree() >= 7]) > 0][0]

    return h.tabulate_lift_evals()


def count_LMFDB_lifts():
    for d in [1,2,3,7,11]:
        evs = load_liftable_from_csv(d)
        print("Total number of forms:", len(evs))
    
def non_steinberg_lift():
    """
    Compute example of non-steinberg lift at p=2 in paper.
    """
    from lmf import db
    cands = [h for h in load_liftable_from_csv(11) if h.p() == 2 and h.level_label() == "23.1"]
    assert len(cands) == 1
    h = cands[0]
    print("Number of lifts:", len(h.get_rational_lifts()))
    lift_reduction_types(h)

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
    print(tabulate(table, headers=['lift', 'q', 'reduction at q', 'Kodaira', 'ord_q(N)', 'semistab'],
                   tablefmt='simple_grid'))
    return table
