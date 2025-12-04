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

class HeckeEig():
    # list of rational lifts of Hecke eigenvalue system
    # can't make it set bc/ lists are not hashable
    def __init__(self, evals, lvl, p):
        self.evals = evals
        self.lvl = lvl
        self.prime = p
        self.rational_lifts = []         

    def __repr__(self):
        return f"Mod {self.p()} Hecke eigenvalue system: {self.eigenvalues()}"

    # Helper methods to make code below cleaner.
    # Not strictly necessary, but prevents us from accidentally
    # changing variables in HeckeEig object.
    def level(self):
        return self.lvl

    def level_label(self):
        return ideal_label(self.level())

    
    def field(self):
        return self.level().number_field()

    def p(self):
        return self.prime
    
    def eigenvalues(self):
        return self.evals

    def add_rational_lift(self, lift):
        self.rational_lifts.append(lift)

    def has_rational_lift(self):
        return self.rational_lifts != []

    def get_rational_lifts(self):
        return self.rational_lifts

    def csv_dict(self):
        if self.has_rational_lift():
            return {'disc': self.field().discriminant(),
                    'evals': self.eigenvalues(),
                    'p': self.p(),
                    'lift_labels': [x[0] for x in self.get_rational_lifts()],
                    'lift_evs': [x[1] for x in self.get_rational_lifts()],
                    }
        else:
            return {'disc': self.field().discriminant(),
                    'evals': self.eigenvalues(),
                    'p': self.p(),
                    }
        
        
        

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


# I think we need to assume that the list of eigenvalues is aligned with the list in lmfdb,
# i.e. that the i-th entry corresponds to the same ideal in both lists. Not happy.

def find_congruent_forms(d, input_file = None, find_all_lifts = False):
    """
    Search LMFDB for Bianchi eigenforms on field on discriminant d
    with eigenvalues congruent to eigenvalue tuples stored in `../data/nonEis_d*.csv`.
    Returns list of hecke eigenvalue objects, whose field `rational_lifts` is populated
    with rational lifts, in the form of a tuple (label, evs, level)
    """

    R.<x> = PolynomialRing(ZZ)
    K.<t> = NumberField(x^2 + d)
    BMFs = db.bmf_forms
    OK = K.maximal_order()
    label = f"2.0.{K.discriminant().abs()}.1"
    if input_file is None:
        input_file = f"../data/nonEis_d{d}.csv"

    # load mod p BMFs from input_file
    # this is small enough that it doesn't make sense to use an iterator. 

    modp_evs = []
    with open(input_file, newline='') as f:
        csv_reader = csv.reader(f, delimiter=';')

        prime_labels = next(csv_reader)[2:]
            
        for row in csv_reader:
            
            lvl = ideal_from_label(K, row[0])
            p = eval(row[1])
            ev = []
            is_valid = True
            for x in row[2:-1]:
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

    # return modp_evs;

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
                print("Found rational lift for", hecke_ev)
                hecke_ev.add_rational_lift({'label': hit['label'], 'evals': hit_ev, 'level': hit_lvl})
                # test level raising condition
                non_LR_primes = []
                for i, pp in enumerate(prime_list):
                    if pp.divides(hit_lvl) and not pp.divides(p*lvl):
                        modp_diff = (ev[i]^2 - (1 + norm(pp))^2) % p
                        if  modp_diff != 0:
                            non_LR_primes.append([modp_diff, prime_label(pp)])
                if non_LR_primes != []:
                    print(f"Lift for Mod {p} ev system {ev} failed level raising condition")
                if p > 3:
                    for val, label in non_LR_primes:
                        print(f"Difference mod {p} is {val} for coeff of {label}")
            
            # TODO: Make sure magma evals are ordered the same way!
            # TODO: Find out exactly what mismatches we should expect.
            # Currently, accept different eigenvalues away from level and p.


    return modp_evs

def string_to_level(K, str):
    """
    Convert the LMFDB level 'str' into an element of the number field K
    """
    
    w = K.gen()
    strlist = list(str)
    if strlist.count("+") == 1:
        plusind = strlist.index("+")
        onecoeff = "".join(strlist[0:plusind])
        wcoeff = "".join(strlist[plusind+1:len(strlist)-1])
        if onecoeff == "":
            onecoeff = 1
        if wcoeff == "":
            wcoeff = 1
            return int(onecoeff) + int(wcoeff)*w
        elif strlist.count("t") == 1:
            return int("".join(strlist[0:len(strlist)-1]))*w
        else:
            return K(int("".join(strlist[0:len(strlist)-1])))


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

def test_write():
    d = 2
    modp_evs = find_congruent_forms(d)
    with open(f'd{d}_liftable.csv', 'w', newline='') as f:
        fieldnames = ['disc','p','evals','lift_labels', 'lift_evs']
        writer = csv.DictWriter(f, delimiter=';', fieldnames=fieldnames)
        writer.writeheader()
        for h in modp_evs:
            if h.has_rational_lift():
                writer.writerow(h.csv_dict())
    print(f"Wrote liftable eigenvalues and lifts to", f'd{d}_liftable.csv')

    with open(f'd{d}_unliftable.csv', 'w', newline='') as f:
        fieldnames = ['disc','p','evals']
        writer = csv.DictWriter(f, delimiter=';', fieldnames=fieldnames)
        writer.writeheader()
        for h in modp_evs:
            if not h.has_rational_lift():
                writer.writerow(h.csv_dict())
    print(f"Wrote unliftable eigenvalues to", f'd{d}_unliftable.csv')

    

    

# def test_search():
#     R.<x> = PolynomialRing(ZZ)
#     d = 3
#     targetEigs = [ 1, 1, 1, 0, 0, 1, 0, 2, 2, 0, 2 ]
#     p = 3

#     K.<w> = NumberField(x^2+d)

#     badPrimeIdeals = [w, 1-3*w]


#     badPrimeNorms = [norm(u) for u in badPrimeIdeals]
#     levelIdeal = prod([K.ideal(u) for u in badPrimeIdeals])


#     labels = []
#     all_lift_primes = []

#     for i in range(len(ll)):
#         if all([u in ll[i]["level_bad_primes"] for u in badPrimeNorms]):
#             currentEigs = ll[i]["hecke_eigs"][0:len(targetEigs)]
#             I = K.ideal(StringToLevel(ll[i]["level_gen"]))
#             ff = I.factor()
#             II = I.intersection(levelIdeal)
#             badPNum = len(II.factor())
#             diffs = [mod(currentEigs[i]-targetEigs[i],p) for i in range(len(currentEigs))]
#             if diffs.count(0) > len(currentEigs) - badPNum - 2 and II == I: 
#                 #if True:
#                 #ll[i]["label"]
#                 labels.append(ll[i]["label"])
#                 (StringToLevel(ll[i]["level_gen"])/levelIdeal).norm()
#                 fac = (II/levelIdeal).factor()
#                 # we add all the primes appearing in the lifting factor to a list 
#                 lift_primes = [fac[i][0].gens_reduced()[0] for i in range( len(fac) )]
#                 for i in range(len(lift_primes)):
#                         if not lift_primes[i] in all_lift_primes:
#                         all_lift_primes.append(lift_primes[i])
#                         #currentEigs
#                         #diffs
#                         # ""





##

# ECs = db.ec_nfcurves


# ranks = []

# for u in labels:

# 	infoEC = {}
# 	LEC = ECs.search({"label": u + "1"}, projection=['rank'], limit=100, info = infoEC)
# 	ranks.append(LEC[0]["rank"])

