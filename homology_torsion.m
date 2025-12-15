
// given the discriminant of an imaginary quadratic field, returns an FP group 
// isomorphic to PGL_2(ZK) (or PSL_2(ZK) if GL is set to false)
// only works for discriminants -3, -4, -7, -8, -11.
BianchiPresentation:=function(disc : GL:=true)

    if GL then 

	if disc eq -3 then 
	    return Group<t, tw, s, j | s^2 = (t*s)^3 = (t*s*tw*s*tw^-1*s*t^-1)^2 = (t^-1*tw*t*s*tw*s*tw^-1*s*t^-1)^3 = ((s*t*s*tw*s*tw^-1*s*t^-1)^-1*t^-1*tw*(s*t*s*tw*s*tw^-1*s*t^-1)*t^-1) = (s*t*s*tw*s*tw^-1*s*t^-1)^-1*t*(s*t*s*tw*s*tw^-1*s*t^-1)*tw = t*tw*t^-1*tw^-1 = j^-1 * t *j * tw * t^-1 = j^-1*tw*j*t^-1 = j^-1 * s * j *(t*s*tw*s*tw^-1*s*t^-1) = 1 >;
	    elif disc eq -4 then 
			return Group<t, tw, s, j | s^2 = (s*tw^-1*s*tw*s*tw^-1)^2 = (s*(s*tw^-1*s*tw*s*tw^-1))^2 = (t*(s*tw^-1*s*tw*s*tw^-1))^2 = (tw*(s*tw^-1*s*tw*s*tw^-1))^2 = (t*s)^3 = (tw*s*(s*tw^-1*s*tw*s*tw^-1))^3 = t*tw*t^-1*tw^-1 = j^-1*t*j*tw = j^-1*tw*j*t^-1 = j^-1*s*j * (tw^-1*s*tw*s*tw^-1) = 1 >;
	    elif disc eq -7 then 
			return Group<t, tw, s, j | s^2 = (t*s)^3 = t*tw*t^-1*tw^-1 = (tw^-1*s*tw*s*t)^2 = j^-1*t*j*t = j^-1*tw*j*tw = j^-1*s*j*s = 1 >;
	    elif disc eq -8 then 
			return Group<t, tw, s, j | s^2 = (t*s)^3 = t*tw*t^-1*tw^-1 = (s*tw*s*tw^-1)^2 = j^-1*t*j*t = j^-1*tw*j*tw = j^-1*s*j*s = 1 >;
	    elif disc eq -11 then 
			return Group<t, tw, s, j | s^2 = (t*s)^3 = t*tw*t^-1*tw^-1 = (tw^-1*s*tw*s*t)^3 = j^-1*t*j*t = j^-1*tw*j*tw = j^-1*s*j*s = 1 >;
	    
	end if;

    else 

		if disc eq -3 then 
		    return Group<t, tw, s | s^2 = (t*s)^3 = (t*s*tw*s*tw^-1*s*t^-1)^2 = (t^-1*tw*t*s*tw*s*tw^-1*s*t^-1)^3 = ((s*t*s*tw*s*tw^-1*s*t^-1)^-1*t^-1*tw*(s*t*s*tw*s*tw^-1*s*t^-1)*t^-1) = (s*t*s*tw*s*tw^-1*s*t^-1)^-1*t*(s*t*s*tw*s*tw^-1*s*t^-1)*tw = t*tw*t^-1*tw^-1 = 1 >;
		    elif disc eq -4 then 
			return Group<t, tw, s | s^2 = (s*tw^-1*s*tw*s*tw^-1)^2 = (s*(s*tw^-1*s*tw*s*tw^-1))^2 = (t*(s*tw^-1*s*tw*s*tw^-1))^2 = (tw*(s*tw^-1*s*tw*s*tw^-1))^2 = (t*s)^3 = (tw*s*(s*tw^-1*s*tw*s*tw^-1))^3 = t*tw*t^-1*tw^-1 = 1 >;
		    elif disc eq -7 then 
			return Group<t, tw, s | s^2 = (t*s)^3 = t*tw*t^-1*tw^-1 = (tw^-1*s*tw*s*t)^2 = 1 >;
		    elif disc eq -8 then 
			return Group<t, tw, s | s^2 = (t*s)^3 = t*tw*t^-1*tw^-1 = (s*tw*s*tw^-1)^2 = 1 >;
		    elif disc eq -11 then 
			return Group<t, tw, s | s^2 = (t*s)^3 = t*tw*t^-1*tw^-1 = (tw^-1*s*tw*s*t)^3 = 1 >;
		end if;

    end if;

		     end function;


ExtraModPClasses:=function(level : GL := true)
    ZK := Order(level);
    K := NumberField(ZK);

    w := K.1;
    T := Matrix(ZK,2,2,[1,1,0,1]);
    Tw := Matrix(ZK,2,2,[1,w,0,1]);
    S := Matrix(ZK,2,2,[0,-1,1,0]);
    // print Discriminant(ZK);
    G := BianchiPresentation(Discriminant(ZK) : GL := GL);

    PL,r := ProjectiveLine(quo<ZK|level>);
    proj_action := function(M,i)
	t,im := r(PL[i]*M,true,false);
	return Index(PL,im);
		   end function;

    action_perms := [[],[],[]];
	
    // we compute the coset action of T, Tw and S	
    for i in [1..#PL] do
	Append(~action_perms[1], proj_action(T,i));
	Append(~action_perms[2], proj_action(Tw,i));
	Append(~action_perms[3], proj_action(S,i));
    end for;

    // we also compute it for J if we want the GL version 
    if GL then 
	if Discriminant(ZK) eq -3 or Discriminant(ZK) eq -4 then 
	    J := Matrix(ZK,2,2,[w,0,0,1]);
	else
			J := Matrix(ZK,2,2,[-1,0,0,1]);
	end if;
	J_perm := [];
	for i in [1..#PL] do 
	    Append(~J_perm, proj_action(J,i));
	end for;
	Append(~action_perms, J_perm);
    end if;

    Symm := Sym(#PL);

    // then we compute the subgroup and return the AQInvariants (piece of the abelian quotient as a f.g. abelian group)
    f := hom<G -> Symm | action_perms >;
    H := sub< G | f >;

    return Sort(AQInvariants(H));
		  end function;

// Compute levels N and primes p where we expect nonlifting of mod p classes of level N.
function ComputeLevelsAndPrimes(d : lowerBound := 1, upperBound:= 1000)
    assert d in [1,2,3,7,11];
    _<x> := PolynomialRing(Rationals());
    // F := NumberField(x^2+d);
    F := QuadFld(d);
    sigma := Automorphisms(F)[2];

    Ids := IdealsUpTo(upperBound,F);
    Ids_no_conj := [];
    for u in Ids do 
	if not sigma(u) in Ids_no_conj and Norm(u) ge lowerBound then 
	    Append(~Ids_no_conj,u);
	end if;
    end for;
    data:=[];

    for u in Ids_no_conj do 
	_,g:=IsPrincipal(u);
	invs:=ExtraModPClasses(u);
	if invs[#invs] ne 0 then 
	    Append(~data,[Eltseq(g), PrimeFactors(invs[#invs])]);
	end if;
    end for;

    return data;

end function;








