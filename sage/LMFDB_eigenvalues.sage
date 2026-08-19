# Uses lmfdb lite, see https://github.com/roed314/lmfdb-lite
# Make sure to install via sage's own pip:
# sage -pip install -U "lmfdb-lite[pgbinary] @ git+https://github.com/roed314/lmfdb-lite.git"


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
    def __init__(self, evals, lvl, p, label=None, is_irrational=False):
        self.evals = evals
        self.lvl = lvl
        self.prime = p
        self.rational_lifts = []
        self.irrational_lifts = []
        self._label = label
        self.is_irrational = is_irrational

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

    def csv_dict(self):
        if self.has_rational_lift():
            return {'level': self.level_label(),
                    'evals': self.eigenvalues(),
                    'p': self.p(),
                    'lift_labels': [x['label'] for x in self.get_rational_lifts()],
                    # 'lift_evs': [x['evals'] for x in self.get_rational_lifts()],
                    }
        else:
            return {'level': self.level_label(),
                    'evals': self.eigenvalues(),
                    'p': self.p(),}

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
        
        

# def compute_magma_modp_Bianchi():
#     """
#     Compute all the mod p non-Eisenstein Bianchi eigenvalue systems,
#     return list [hecke_eig1, hecke_eig2, ...]
#     """
    
#     if 'computations/sage' not in os.getcwd():
#         print("Ensure that the current directory is BMF/computations/sage, otherwise magma won't find the spec!")
#         return 1

#     # change working directory for magma for readability
#     if "test" not in magma.eval('System("pwd")'):
#         magma.chdir("../../test")
    
#     magma.load("loading_script_modp_nonEis.m")

#     # magma.eval("TestNonEis();")
#     print("Finished computing non-Eisenstein systems, written to ../data/nonEis_d2.csv")
#     return 0

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
        input_file = f"../data/nonEis_d{d}.csv"

    # load mod p BMFs from input_file
    # this is small enough that it doesn't make sense to use an iterator. 

    modp_evs = []
    with open(input_file) as f:
        csv_reader = csv.reader(f, delimiter=';')

        prime_labels = next(csv_reader)[2:-1]
            
        for row in csv_reader:
            if row == []:
                continue
            lvl = ideal_from_label(K, row[0])
            p = eval(row[1])
            ev = []
            is_valid = True
            for i, x in enumerate(row[2:-1]):
                try:
                    ev.append(ZZ(x))
                except TypeError:
                    if x == "-":
                        ev.append(ZZ(0))
                    else:
                        print(f"Failed to coerce {x} to integer")
                        is_valid = False
                        break
            if is_valid:
                modp_evs.append(HeckeEig(ev, lvl, p))


    max_ev_len = max(len(h.eigenvalues()) for h in modp_evs)
    prime_list = [prime_from_label(K,lab) for lab in prime_labels]

    print(f"Successfully loaded {len(modp_evs)} eigenvalues from {input_file}")

    hits = BMFs.search({'dimension': 1, "field_label": label}, ['label', 'hecke_eigs', 'level_bad_primes', 'level_label'], limit = SEARCH_COUNT)
    print("Successfully pulled Hecke eigenvalue systems from LMFDB")

    for hit in hits:
        hit_lvl = ideal_from_label(K, hit['level_label'])
        hit_ev = hit['hecke_eigs']
        

        for hecke_ev in modp_evs:
            # unless we're trying to find all lifts, skip to next once we have a rational lift
            if not find_all_lifts and hecke_ev.has_rational_lift():
                continue

            lvl = hecke_ev.level()
            # if the level of potential lift isn't divisible by level of hecke_ev, skip to next
            if not lvl.divides(hit_lvl):
                continue
            
            ev = hecke_ev.eigenvalues()
            p = hecke_ev.p()


            min_len = min(len(ev), len(hit_ev))
            # We allow a lift to have different eigenvalues corresponding to "bad" primes
            # This includes:
            #  - Primes dividing p,
            #  - Primes above the level of the lift
            #  - Primes dividing the discriminant of the field

            if all(((ev[i] - hit_ev[i]) % p == 0 or prime_list[i].divides(p*hit_lvl))
                   for i in range(min_len)):
                print("Found rational lift for", hecke_ev, "of level", ideal_label(hit_lvl))
                hecke_ev.add_rational_lift(hit['label'])
            
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

def test_write(d, modp_evs = None):
    if modp_evs is None: 
        modp_evs = find_congruent_forms(d, find_all_lifts=true)

    with open(f'd{d}_liftable.csv', 'w', newline='') as f:
        fieldnames = ['level','p','evals','lift_labels']
        writer = csv.DictWriter(f, delimiter=';', fieldnames=fieldnames)
        writer.writeheader()
        for h in modp_evs:
            if h.has_rational_lift():
                writer.writerow(h.csv_dict())
    print(f"Wrote liftable eigenvalues and lifts to", f'd{d}_liftable.csv')

    with open(f'd{d}_unliftable.csv', 'w', newline='') as f:
        fieldnames = ['level','p','evals']
        writer = csv.DictWriter(f, delimiter=';', fieldnames=fieldnames)
        writer.writeheader()
        for h in modp_evs:
            if not h.has_rational_lift():
                writer.writerow(h.csv_dict())
    print(f"Wrote unliftable eigenvalues to", f'd{d}_unliftable.csv')

def load_unliftable_from_csv(d):
    K = QuadFld(d)
    modp_evs = []
    with open(f'd{d}_unliftable.csv', 'r', newline='') as f:
        # fieldnames = ['level','p','evals','lift_labels', 'lift_evs']
        reader = csv.DictReader(f, delimiter=';')
        for r in reader:
            lvl = ideal_from_label(K, r['level'])
            modp_evs.append(HeckeEig(eval(r['evals']), lvl, ZZ(r['p'])))
    return modp_evs

def load_liftable_from_csv(d):
    K = QuadFld(d)
    modp_evs = []
    with open(f'd{d}_liftable.csv', 'r', newline='') as f:
        # fieldnames = ['level','p','evals','lift_labels', 'lift_evs']
        reader = csv.DictReader(f, delimiter=';')
        for r in reader:
            lvl = ideal_from_label(K, r['level'])
            h = HeckeEig(eval(r['evals']), lvl, ZZ(r['p']))
            for lift in eval(r['lift_labels']):
                h.add_rational_lift(lift)
            modp_evs.append(h)
            
    return modp_evs

def load_irrat_results(d):
    """Load list of all irrational forms with lifts from ../data/irrat_lifts_d{d}.csv,
    ignoring the ones where we haven't found anything yet."""
    K = QuadFld(d)
    irrat_lifts = []
    with open(f'../data/irrat_lifts_d{d}.csv') as f:
        reader = csv.DictReader(f, delimiter=';')
        for r in reader:
            if not "No" in r[' lift level']:
                p = ZZ(r[' prime'])
                lvl = ideal_from_label(K, r['level'])
                h = HeckeEig(None, lvl, p)
                h.add_irrational_lift(r[' lift level'])
                irrat_lifts.append(h)

    return irrat_lifts
        

def load_all_from_csv(d):
    """
    Load liftable and unliftable forms, including those where we haven't found lifts yet.
    """
    # first, collect forms from three different sources
    rat_l_forms = load_liftable_from_csv(d)
    rat_ul_forms = load_unliftable_from_csv(d)
    irrat_forms = load_irrat_forms(d)

    # then eliminate the ones where we found a lift using magma
    irrat_lifts = load_irrat_lifts(d)

    for h in rat_ul_forms:
        cand = [g for g in irrat_lifts if h.level() == g.level() and h.p() == g.p()]
        if cand:
            h.add_irrational_lift(cand[0].get_irrational_lifts()[0])

    for h in irrat_forms:
        cand = [g for g in irrat_lifts if h.level() == g.level() and h.p() == g.p()]
        if cand:
            h.add_irrational_lift(cand[0].get_irrational_lifts()[0])

    return rat_l_forms + rat_ul_forms + irrat_forms


def load_irrat_forms(d):
    """
    Load irrational mod p Hecke eigenvalue systems 
    """
    K = QuadFld(d)
    input_file = f"../data/nonEis_d{d}.csv"

    irrat_forms = []
    with open(input_file, newline='') as f:
        csv_reader = csv.reader(f, delimiter=';')
        # TODO: find a good way to parse Hecke evals
        for i, row in enumerate(csv_reader):
            if i == 0:
                continue
            p = eval(row[1])
            
            lvl = ideal_from_label(K, row[0])
            if any("[" in r for r in row):
                irrat_forms.append(HeckeEig(
                    None, lvl, p, is_irrational=True
                ))

    return irrat_forms
    
def ul_filter_csv(d):
    """
    Create magma-readable file with EBMFs which don't lift in the lmfdb
    """
    F = QuadFld(d)
    input_file = f"../data/nonEis_d{d}.csv"
    output_file = f"../data/lmfdbNonlift_d{d}.csv"
        
    modp_evs = load_unliftable_from_csv(d)
    modp_evs = [h for h in modp_evs]
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
            for x in row[2:-1]:
                try:
                    evs.append(ZZ(0) if x == "-" else ZZ(x))
                except TypeError:
                    print(f"Failed to coerce {x} to integer in row {row[0]};{row[1]}")
                    is_valid = False
                    break
            if not is_valid:
                continue
            if any(row[0] == h.level_label() and p == h.p() and list(h.eigenvalues()) == evs
                   for h in modp_evs):
                lines.append(";".join(row))
        
    with open(output_file, "w") as f:
        f.write("\n".join(lines))
    

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
        
def count_l_vs_ul(d, norm_bd=500, p_bd=20, print_unliftable=False):
    all_evs = load_all_from_csv(d)
    ls = [h for h in all_evs if h.has_lift()]
    uls = [h for h in all_evs if not h.has_lift()]
    ls = [ev for ev in ls if ev.level().norm() <= norm_bd and ev.p() < p_bd]
    uls = [ev for ev in uls if ev.level().norm() <= norm_bd and ev.p() < p_bd]
    print(f"{len(ls)}/{len(ls) + len(uls)} of forms with p < {p_bd} and lvl norm < {norm_bd} have char 0 lifts")
    if print_unliftable and len(uls) > 0:
        print("Forms without known lift:" + " "*(len(f"{uls[0]}") - 18) + "Irrational?")
        for ev in uls:
            print(f"{ev}\t {ev.is_irrational}")
    return len(ls)/(len(ls) + len(uls))

    
def count_single_prime_lifts(d):
    F = QuadFld(d)
    ls = load_liftable_from_csv(d)
    def is_good(h):
        newlvl = level_from_BMF_label(F, h.get_rational_lifts()[0])/h.level()
        return newlvl.is_prime() and newlvl.is_coprime(h.p())
    return len([h for h in ls if is_good(h)])/len(ls)
