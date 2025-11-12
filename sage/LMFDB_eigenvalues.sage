# Uses lmfdb lite, see https://github.com/roed314/lmfdb-lite
# Make sure to install via sage's own pip:
# sage -pip install -U "lmfdb-lite[pgbinary] @ git+https://github.com/roed314/lmfdb-lite.git"


from lmf import db
import re

DEBUG_LEVEL = ZZ(1e5)

# NOTE: We cannot search on Hecke eigenvalue modulo p because the 'hecke_eigs' field is jsonb, not array.
# So we have to search through all the forms. 

# TODO: make this exportable to csv
# TODO: make this loadable from csv
class hecke_eig():
    evals = []
    has_rational_lift = False
    rational_lifts = []
    def __init__(self, evals, field, ):
        self.evals = evals
        self.field = field


x = hecke_eig([0,1,2], K)
x.has_rational_lift

def compute_magma_modp_Bianchi():
    """
    Compute all the mod p non-Eisenstein Bianchi eigenvalue systems,
    return list [hecke_eig1, hecke_eig2, ...]
    """
    magma.attach_spec("../../spec")
    magma.load


# I think we need to assume that the list of eigenvalues is aligned with the list in lmfdb,
# i.e. that the i-th entry corresponds to the same ideal in both lists. Not happy.

def find_congruent_form(K, eval_list, p):
    """
    Search LMFDB for Bianchi eigenforms of level `lvl` on field `K` with eigenvalues congruent to
    eigenvalue tuples in `eval_list` (of the form [[f1_evals, f1_level], ...]) mod `p`.
    Returns hits in the form of a dictionary { 'eval' : [('label', 'hecke_eigs', 'level')]}
    """
    
    BMFs = db.bmf_forms
    OK = K.maximal_order()
    a = K.gen()                 # necessary for parsing
    label = f"2.0.{K.discriminant().abs()}.1"

    hits = BMFs.search({'dimension': 1, "field_label": label}, ['label', 'hecke_eigs', 'level_bad_primes', 'level_ideal'], limit = DEBUG_LEVEL)

    output_dict = {}
    for hit in hits:
        # add in multiplications in ideal strings
        hit_lvl = K(re.sub(r"(\d+)([a-zA-Z])",
                           r"\1*\2", hit['level_ideal']))
        hit_ev = hit['hecke_eigs']

        for ev, lvl in eval_list:
            # print([(ev[i] - hit_ev[i]) %p for i in range(min(len(ev),len(hit_ev))) ])
            # print(lvl.divides(hit_lvl))
            
            mismatches = [(i,ev[i],hit_ev[i])
                          for i in range(min(len(ev),len(hit_ev)))
                          if ev[i] % p != hit_ev[i] % p]
            # TODO: find out exactly how many mismatches we should have
            # currently, away from level and characteristic
            if lvl.divides(hit_lvl) and len(mismatches) <= len(factor(hit_lvl))+2:
                try:
                    # consider if there's a better way to index our forms?
                    # lists are not hashable but tuples are
                    output_dict[str(ev)].append((hit['label'], str(hit_ev), hit_lvl))
                except KeyError:
                    output_dict[str(ev)] = [(hit['label'], str(hit_ev), hit_lvl)]
                    
    return output_dict
    
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


def test():
    # TODO: make sure we use same generators as LMFDB!
    R.<x> = PolynomialRing(ZZ)
    d = 3
    p = 3
    K.<w> = NumberField(x^2+d)
    lvl =  K.ideal(w, (1-3*w))
    eval_list = [([ 1, 1, 1, 0, 0, 1, 0, 2, 2, 0, 2 ], lvl)]
    return find_congruent_form(K, eval_list, p)

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

