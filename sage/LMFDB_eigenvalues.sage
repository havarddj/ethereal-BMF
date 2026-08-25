# Uses lmfdb lite, see https://github.com/roed314/lmfdb-lite
# Make sure to install via sage's own pip:
# sage -pip install -U "lmfdb-lite[pgbinary] @ git+https://github.com/roed314/lmfdb-lite.git"
# and tabulate: sage -pip install -U tabulate

from lmf import db
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
    def __init__(self, evals, lvl, p, label=None, eig_dim = None, gen_dim = None):
        self.evals = evals
        self.lvl = lvl
        self.prime = p
        self.rational_lifts = []
        self.irrational_lifts = []
        self._label = label
        self.eig_dim = eig_dim
        self.gen_dim = gen_dim


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

    def get_irrational_lifts(self):
        return self.irrational_lifts

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

    R.<x> = PolynomialRing(ZZ)
    K = QuadFld(d)
    BMFs = db.bmf_forms
    OK = K.maximal_order()
    label = f"2.0.{K.discriminant().abs()}.1"
    if input_file is None:
        input_file = f"../data/nonEis_d{d}_v2.csv"

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
    """ Load all Hecke eigenvalue systems from nonEis_dX_v2.csv files
    """
    K = QuadFld(d)
    
    input_file = f"../data/nonEis_d{d}_v2.csv"
    
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
    irrat_file= f'../data/irrat_lifts_d{d}_v2.csv'
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
                h = HeckeEig(evals, lift_lvl, 0)
                irrat_lifts.append(h)
    # Since we don't store the exact one, we need to match up with nonliftable:
    if not modp_evs:
        modp_evs = load_unliftable_from_csv(d)

    undetected = irrat_lifts
    for g in modp_evs:
        cands = [h for h in irrat_lifts if g.level().divides(h.level()) and is_irrational_lift(h, g)]
        if len(cands) > 1:
            print("WARNING: Multiple candidate lifts available for", g)
            print("Levels:", [c.level_label() for c in cands])
        for c in cands:
            g.add_irrational_lift(c)
            undetected.remove(c)
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
    Load liftable and unliftable forms, including those where we haven't found lifts yet.
    """
    # first, collect forms from three different sources
    lifties = load_liftable_from_csv(d)
    unlifties = load_unliftable_from_csv(d)
    _ = load_irrat_lifts(d, modp_evs = unlifties)
    
    return lifties + unlifties

    
def ul_filter_csv(d):
    """
    Create magma-readable file with EBMFs which don't lift in the lmfdb
    """
    F = QuadFld(d)
    input_file = f"../data/nonEis_d{d}_v2.csv"
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
    from tabulate import tabulate
    def make_table(lifts):
        return [["Ethereal forms total:", len(lifts)],
                    ["Ethereal forms with lift:", len([h for h in lifts if h.has_lift()])],
                    ["Ethereal Forms with no lift:", len([h for h in lifts if not h.has_lift()])],
                    ["Fp-valued forms total:", len([h for h in lifts if h.is_rational()])],
                    ["Fp-valued forms with rational lift:", len([h for h in lifts if h.is_rational() and h.has_rational_lift()])],
                    ["Fp-valued forms with irrational lift:", len([h for h in lifts if h.is_rational() and h.has_irrational_lift()])],
                    ["Fp-valued forms with lift (total):", len([h for h in lifts if h.is_rational() and h.has_lift()])],
                    ["Fp-valued forms with no lift:", len([h for h in lifts if not h.has_lift() and h.is_rational()])],
                ]
    total_lifts = []
    total_small_lifts = []
    total_small_p_lifts = []
    for d in [1,2,3,7,11]:
        all_lifts = load_all_from_csv(d)
        small_p_lifts = [h for h in all_lifts if 2 < h.p() < 20]
        small_lifts = [h for h in all_lifts if 2 < h.p() < 20 and h.level().norm() <= 500]
        total_lifts += all_lifts
        total_small_p_lifts += small_p_lifts
        total_small_lifts += small_lifts
        print("\n")
        print(tabulate(make_table(all_lifts), headers=[f"d={d}, no conditions", "Count"], tablefmt="grid"), "\n")
        print(tabulate(make_table(small_p_lifts), headers=[f"d={d}, 2 < p < 20", "Count"], tablefmt="simple_grid"))
        print(tabulate(make_table(small_lifts), headers=[f"d={d}, 2 < p < 20, Nm(n) <= 500", "Count"], tablefmt="simple_grid"))
        print("\n")
    print("\n")
    print(tabulate(make_table(total_lifts), headers=[f"All d, no conditions", "Count"], tablefmt="double_grid"), "\n")
    print(tabulate(make_table(total_small_p_lifts), headers=[f"All d, 2 < p < 20", "Count"], tablefmt="double_grid"))
    print(tabulate(make_table(total_small_lifts), headers=[f"All d, 2 < p < 20, Nm(n) <= 500", "Count"], tablefmt="double_grid"))
    print("\n")
    
