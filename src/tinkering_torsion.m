
AttachSpec("../BMF/spec");
load "homology_torsion.m";

// given the discriminant of an imaginary quadratic field, returns an FP group 
// isomorphic to PGL_2(ZK) (or PSL_2(ZK) if GL is set to false)
// only works for discriminants -3, -4, -7, -8, -11.
BianchiPresentationNew:=function(disc : GL:=true)

    if GL then 

	if disc eq -3 then 
	    return Group<t, tw, s, j | j^2 * (s*t*s*tw*s*tw^-1*s*t^-1)^-1 = s^2 = (t*s)^3 = (t*s*tw*s*tw^-1*s*t^-1)^2 = (t^-1*tw*t*s*tw*s*tw^-1*s*t^-1)^3 = ((s*t*s*tw*s*tw^-1*s*t^-1)^-1*t^-1*tw*(s*t*s*tw*s*tw^-1*s*t^-1)*t^-1) = (s*t*s*tw*s*tw^-1*s*t^-1)^-1*t*(s*t*s*tw*s*tw^-1*s*t^-1)*tw = t*tw*t^-1*tw^-1 = j^-1 * t *j * tw * t^-1 = j^-1*tw*j*t^-1 = j^-1 * s * j *(t*s*tw*s*tw^-1*s*t^-1) = 1 >;
    elif disc eq -4 then 
		return Group<t, tw, s, j | j^2 * (s*tw^-1*s*tw*s*tw^-1)^-1 = s^2 = (s*tw^-1*s*tw*s*tw^-1)^2 = (s*(s*tw^-1*s*tw*s*tw^-1))^2 = (t*(s*tw^-1*s*tw*s*tw^-1))^2 = (tw*(s*tw^-1*s*tw*s*tw^-1))^2 = (t*s)^3 = (tw*s*(s*tw^-1*s*tw*s*tw^-1))^3 = t*tw*t^-1*tw^-1 = j^-1*t*j*tw = j^-1*tw*j*t^-1 = j^-1*s*j * (tw^-1*s*tw*s*tw^-1) = 1 >;
    elif disc eq -7 then 
		return Group<t, tw, s, j | j^2 = s^2 = (t*s)^3 = t*tw*t^-1*tw^-1 = (tw^-1*s*tw*s*t)^2 = j^-1*t*j*t = j^-1*tw*j*tw = j^-1*s*j*s = 1 >;
    elif disc eq -8 then 
		return Group<t, tw, s, j | j^2 = s^2 = (t*s)^3 = t*tw*t^-1*tw^-1 = (s*tw*s*tw^-1)^2 = j^-1*t*j*t = j^-1*tw*j*tw = j^-1*s*j*s = 1 >;
    elif disc eq -11 then 
		return Group<t, tw, s, j | j^2 = s^2 = (t*s)^3 = t*tw*t^-1*tw^-1 = (tw^-1*s*tw*s*t)^3 = j^-1*t*j*t = j^-1*tw*j*tw = j^-1*s*j*s = 1 >;
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


ExtraModPClassesNew:=function(level : GL := true)
    ZK := Order(level);
    K := NumberField(ZK);

    w := K.1;
    T := Matrix(ZK,2,2,[1,1,0,1]);
    Tw := Matrix(ZK,2,2,[1,w,0,1]);
    S := Matrix(ZK,2,2,[0,-1,1,0]);
    // print Discriminant(ZK);
    G := BianchiPresentationNew(Discriminant(ZK) : GL := GL);

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


procedure CompareOldAndNew(d,upperBd)
	K := QuadFld(d);
	for I in IdealsUpTo(upperBd,K) do 
		A_old := ExtraModPClasses(I);
		A_new := ExtraModPClassesNew(I);
		m_old := Max(A_old);
		m_new := Max(A_new);

		if m_old ne 0 and m_new ne 0 then 
			p_old := PrimeFactors(m_old);
			p_new := PrimeFactors(m_new);
			missing := [p : p in p_new | not p in p_old];
			if #missing ne 0 then 
				print LMFDBLabel(I), missing;
			end if;
		end if;
	end for;
end procedure;


for d in [1,2,3,7,11] do 
	printf "Comparing old and new for d = %o\n", d;
	print "The following levels missed the following primes:";
	CompareOldAndNew(d,1000);
	print "";
end for;





///////////


AttachSpec("~/Dropbox/Bianchi mod p periods/CODE/!archive/KleinianGroups/klngpspec");

function QuatToMatrix(g)
	K := BaseField(Parent(g));
	return MatrixRing(K,2) ! [g[1]+g[2],g[3]+g[4],g[3]-g[4],g[1]-g[2]];
end function;

function MatrixToQuat(m,O)
	B<i,j,k> := Algebra(O);
	return (m[1,1]*(1+i)+m[2,2]*(1-i)+m[1,2]*(j+k)+m[2,1]*(j-k))/2;
end function;



for d in [1,2,3,7,11] do 

	
	O := BianchiOrder(-d);
	_, F, E := NormalizedBasis(O : GroupType := "Units");
	P, gens := Presentation(F, E, O);   // P ≅ PGL2(O_K), gens are elements of O

	PP := Simplify(P);
	G := BianchiPresentationNew(Discriminant(CoefficientRing(O)));

	effort := 1;
	ran_out := false;
	is_isom := false;
	while not is_isom do 
		is_isom := SearchForIsomorphism(PP,G,effort);
		effort +:= 1;
		if effort eq 15 then 
			ran_out := true;
			break;
		end if;
	end while;

	if ran_out then
		printf "Ran out of effort trying to find isomorphism for d = %o\n";
	else 
		if is_isom then 
			printf "Found isomorphism for d = %o with effort %o\n", d, effort;
		else
			printf "Couldn't find isomorphism for d = %o", d;
		end if;
	end if;
end for;





